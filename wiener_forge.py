#!/usr/bin/env python3
"""
Wiener's Attack on RSA — outil pédagogique de cryptanalyse
===========================================================
Démonstration de l'attaque de Wiener (1990) :
si d < (1/3) * N^(1/4), alors d peut être retrouvé
par développement en fractions continues de e/N.

Auteur : Rachid Ghodbane - rvsh0x

Usage :
    python3 wiener_rsa.py <e_hex_ou_decimal> <N_hex_ou_decimal>

Exemple :
    python3 wiener_rsa.py 0xf70b3b...  0x0207a7...

Dépendances :
    pip install pycryptodome

Note : L'attaque de Wiener nécessite SageMath pour la partie
       fractions continues. Ce script inclut une implémentation
       pure Python de l'algorithme.
"""

import sys
import math
from Crypto.PublicKey import RSA


# ─────────────────────────────────────────────────────────────
# 1. FRACTIONS CONTINUES (implémentation pure Python)
# ─────────────────────────────────────────────────────────────

def continued_fraction(numerator, denominator):
    """
    Calcule le développement en fraction continue de numerator/denominator
    via l'algorithme d'Euclide.

    Retourne la liste des quotients partiels [a0, a1, a2, ...]
    tels que :
        numerator/denominator = a0 + 1/(a1 + 1/(a2 + ...))
    """
    cf = []
    while denominator:
        q = numerator // denominator
        cf.append(q)
        numerator, denominator = denominator, numerator - q * denominator
    return cf


def convergents(cf):
    """
    Calcule toutes les convergentes (approximations successives) à partir
    d'un développement en fraction continue.

    Pour [a0; a1, a2, a3, ...] on obtient les fractions :
        h_n / k_n  où :
            h_n = a_n * h_{n-1} + h_{n-2}
            k_n = a_n * k_{n-1} + k_{n-2}

    Retourne une liste de couples (numérateur, dénominateur).
    """
    n_prev, n_curr = 1, cf[0]     # h_{-1} = 1, h_0 = a_0
    d_prev, d_curr = 0, 1         # k_{-1} = 0, k_0 = 1

    yield (n_curr, d_curr)

    for a in cf[1:]:
        n_prev, n_curr = n_curr, a * n_curr + n_prev
        d_prev, d_curr = d_curr, a * d_curr + d_prev
        yield (n_curr, d_curr)


# ─────────────────────────────────────────────────────────────
# 2. ATTAQUE DE WIENER
# ─────────────────────────────────────────────────────────────

def wiener_attack(e, N):
    """
    Tente de retrouver d par l'attaque de Wiener.

    Principe :
        - On sait que e*d = 1 + k*phi(N), donc e/N ≈ k/d
        - k/d est une convergente du développement en fraction
          continue de e/N si d < (1/3) * N^(1/4)
        - On teste chaque dénominateur comme candidat pour d

    Retourne d si trouvé, None sinon.
    """
    print("[*] Calcul du développement en fractions continues de e/N...")
    cf = continued_fraction(e, N)
    print(f"    {len(cf)} quotients partiels calculés.")

    # Message test pour valider chaque candidat d
    m_test = 12345
    c_test = pow(m_test, e, N)

    print("[*] Parcours des convergentes...")
    for i, (k, d_candidate) in enumerate(convergents(cf)):
        if d_candidate == 0:
            continue

        # Test : si d_candidate est le vrai d, il déchiffre c_test en m_test
        try:
            m_recovered = pow(c_test, d_candidate, N)
            if m_recovered == m_test:
                print(f"[+] d trouvé à la convergente #{i} !")
                print(f"    d = {d_candidate}")
                return d_candidate
        except Exception:
            continue

    print("[-] Attaque de Wiener échouée : d est probablement trop grand.")
    print("    Vérifiez que d < (1/3) * N^(1/4).")
    return None


# ─────────────────────────────────────────────────────────────
# 3. FACTORISATION DE N À PARTIR DE d
# ─────────────────────────────────────────────────────────────

def isqrt(n):
    """Racine carrée entière (Newton)."""
    if n < 0:
        raise ValueError("Racine carrée d'un nombre négatif.")
    if n == 0:
        return 0
    x = n
    y = (x + 1) // 2
    while y < x:
        x, y = y, (y + n // y) // 2
    return x


def factor_n_from_d(e, d, N):
    """
    Retrouve p et q à partir de (e, d, N).

    Algorithme :
        1. On sait que e*d - 1 = k * phi(N)
        2. phi(N) = N - p - q + 1, donc p+q = N - phi(N) + 1
        3. p et q sont les racines de X² - (p+q)X + N = 0
        4. On teste les diviseurs de e*d - 1 pour trouver phi(N)
           tel que le discriminant soit un carré parfait.

    Retourne (p, q) ou (None, None) si échec.
    """
    print("[*] Factorisation de N à partir de d...")

    k_phi = e * d - 1  # = k * phi(N)

    # On cherche un diviseur de k_phi qui soit un phi(N) valide
    # Stratégie : tester les petits diviseurs de k_phi
    # (en pratique k est petit, souvent < quelques milliers)

    # Méthode directe : tester k de 1 à une borne raisonnable
    # phi(N) = k_phi / k doit vérifier :
    #   - p+q = N - phi(N) + 1 est positif et < N
    #   - discriminant (p+q)^2 - 4N est un carré parfait

    for k in range(1, 1_000_000):
        if k_phi % k != 0:
            continue

        phi_candidate = k_phi // k

        # p+q doit être entier et raisonnable
        s = N - phi_candidate + 1  # p + q
        if s <= 0:
            continue

        # Discriminant de X² - s*X + N = 0
        discriminant = s * s - 4 * N
        if discriminant < 0:
            continue

        # Vérifier que le discriminant est un carré parfait
        sqrt_disc = isqrt(discriminant)
        if sqrt_disc * sqrt_disc != discriminant:
            continue

        p = (s + sqrt_disc) // 2
        q = (s - sqrt_disc) // 2

        if p * q == N and p > 1 and q > 1:
            print(f"[+] Facteurs trouvés pour k = {k} !")
            print(f"    p = {p}")
            print(f"    q = {q}")
            return p, q

    print("[-] Factorisation directe échouée, tentative probabiliste...")
    return factor_probabilistic(e, d, N)


def factor_probabilistic(e, d, N):
    """
    Algorithme probabiliste de Miller-Rabin pour factoriser N
    à partir de (e, d).

    Repose sur le fait que e*d - 1 = 2^s * t avec t impair,
    et que pour g aléatoire, g^(2^i * t) mod N peut révéler
    un facteur de N via le PGCD.
    """
    import random

    k_phi = e * d - 1

    # Écrire k_phi = 2^s * t avec t impair
    s = 0
    t = k_phi
    while t % 2 == 0:
        t //= 2
        s += 1

    print(f"    k*phi(N) = 2^{s} × t (t impair)")
    print("    Recherche probabiliste en cours...")

    for attempt in range(100):
        g = random.randint(2, N - 1)
        x = pow(g, t, N)

        if x == 1 or x == N - 1:
            continue

        for _ in range(s - 1):
            y = pow(x, 2, N)
            if y == 1:
                p = math.gcd(x - 1, N)
                if 1 < p < N:
                    q = N // p
                    if p * q == N:
                        print(f"[+] Facteurs trouvés (tentative #{attempt+1}) !")
                        print(f"    p = {p}")
                        print(f"    q = {q}")
                        return p, q
            x = y

        p = math.gcd(x - 1, N)
        if 1 < p < N:
            q = N // p
            if p * q == N:
                print(f"[+] Facteurs trouvés (tentative #{attempt+1}) !")
                print(f"    p = {p}")
                print(f"    q = {q}")
                return p, q

    print("[-] Factorisation probabiliste échouée.")
    return None, None


# ─────────────────────────────────────────────────────────────
# 4. CONSTRUCTION DE LA CLÉ PRIVÉE RSA
# ─────────────────────────────────────────────────────────────

def build_private_key(N, e, d, p, q, output_file="id_rsa"):
    """
    Construit et sauvegarde la clé privée RSA au format PEM.

    Une clé PKCS#1 contient :
        - N, e, d      : paramètres RSA de base
        - p, q         : facteurs premiers de N
        - dp = d mod (p-1)  )
        - dq = d mod (q-1)  ) composantes CRT pour accélérer le déchiffrement
        - qinv = q^-1 mod p )

    pycryptodome calcule les composantes CRT automatiquement.
    """
    print(f"[*] Construction de la clé privée RSA...")

    try:
        rsa_key = RSA.construct((N, e, d, p, q))
        pem = rsa_key.export_key("PEM")

        with open(output_file, "wb") as f:
            f.write(pem)

        print(f"[+] Clé privée sauvegardée : {output_file}")
        print(f"    Taille de la clé : {rsa_key.size_in_bits()} bits")
        print()
        print("    Pour l'utiliser en SSH :")
        print(f"      chmod 600 {output_file}")
        print(f"      ssh -i {output_file} user@host")

        return pem

    except Exception as ex:
        print(f"[-] Erreur lors de la construction de la clé : {ex}")
        return None


# ─────────────────────────────────────────────────────────────
# 5. PIPELINE PRINCIPAL
# ─────────────────────────────────────────────────────────────

def parse_int(s):
    """Accepte un entier en décimal ou en hexadécimal (préfixe 0x)."""
    s = s.strip()
    if s.lower().startswith("0x"):
        return int(s, 16)
    return int(s)


def attack(e, N, output_file="id_rsa"):
    """
    Pipeline complet :
        1. Attaque de Wiener → d
        2. Factorisation de N → p, q
        3. Construction de la clé privée → id_rsa
    """
    print("=" * 60)
    print("  Attaque de Wiener sur RSA")
    print("=" * 60)
    print(f"  N = {N}")
    print(f"  e = {e}")
    print(f"  N ({N.bit_length()} bits), e ({e.bit_length()} bits)")
    print()

    # Vérification préliminaire
    bound = (N ** (1/4)) / 3
    print(f"[*] Borne de Wiener : d doit être < {bound:.2e} pour que l'attaque fonctionne.")
    print()

    # Étape 1 : Retrouver d
    d = wiener_attack(e, N)
    if d is None:
        print("\n[!] Abandon : d introuvable par fractions continues.")
        return None
    print()

    # Étape 2 : Factoriser N
    p, q = factor_n_from_d(e, d, N)
    if p is None:
        print("\n[!] Abandon : impossible de factoriser N.")
        return None
    print()

    # Vérifications de cohérence
    assert p * q == N, "Erreur : p * q ≠ N"
    assert pow(pow(42, e, N), d, N) == 42, "Erreur : e et d ne sont pas inverses"
    print("[+] Vérifications de cohérence : OK")
    print()

    # Étape 3 : Construire la clé privée
    pem = build_private_key(N, e, d, p, q, output_file)
    print()
    print("=" * 60)
    print("  Attaque terminée avec succès.")
    print("=" * 60)
    return pem


# ─────────────────────────────────────────────────────────────
# 6. POINT D'ENTRÉE
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        print("\nUsage : python3 wiener_rsa.py <e> <N> [output_file]")
        print("  e et N peuvent être en décimal ou en hexadécimal (0x...)")
        sys.exit(1)

    e_arg = parse_int(sys.argv[1])
    N_arg = parse_int(sys.argv[2])
    out   = sys.argv[3] if len(sys.argv) > 3 else "id_rsa"

    attack(e_arg, N_arg, out)
