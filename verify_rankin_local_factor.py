#!/usr/bin/env python3
"""Exact Rankin local-factor and character-direction algebra.

Standard library only. Polynomials have integer coefficients in A,B,C,D,X.
The rational generating-function identity is proved by cross multiplication
of the Binet denominators (A-B)(C-D), rather than finite coefficient tests.
Character identities are checked in the free character-exponent lattice.
The ordered Mok kernel and its analytic endpoint normalization remain
mathematical inputs: this verifier does not prove either identification.
"""

import argparse
import json
from pathlib import Path
import sys


VARIABLES = ("A", "B", "C", "D", "X")
ZERO = (0, 0, 0, 0, 0)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def add(left, right):
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, 0) + coefficient
    return {m: c for m, c in result.items() if c}


def scale(poly, scalar):
    return {m: scalar * c for m, c in poly.items() if scalar * c}


def multiply(left, right):
    result = {}
    for m, c in left.items():
        for n, d in right.items():
            monomial = tuple(a + b for a, b in zip(m, n))
            result[monomial] = result.get(monomial, 0) + c * d
    return {m: c for m, c in result.items() if c}


def variable(index):
    return {tuple(int(j == index) for j in range(5)): 1}


def serialize(poly):
    return [
        {"exponents": list(monomial), "coefficient": coefficient}
        for monomial, coefficient in sorted(poly.items())
    ]


def verify_character_directions():
    # Coordinates are eta, psi, omega_h. Addition means multiplication
    # of characters, and inversion negates all coordinates. Checking in
    # the free lattice implies the identities for every compatible
    # specialization to finite-order characters.
    def product(*characters):
        return tuple(sum(entries) for entries in zip(*characters))

    def power(character, exponent):
        return tuple(exponent * entry for entry in character)

    eta = (1, 0, 0)
    psi = (0, 1, 0)
    omega_h = (0, 0, 1)
    chi = product(power(eta, -1), psi)
    omega_f = product(eta, psi, omega_h)
    omega_f_eta = product(omega_f, power(eta, -2))
    require(omega_f_eta == product(chi, omega_h),
            "The eta^-1 target-twist central character is incorrect.")
    # Conjugating a finite-order central character inverts it. Thus the
    # Rankin pair f_eta^rho,h has determinant character chi^-1.
    rankin_character = product(power(omega_f_eta, -1), omega_h)
    require(rankin_character == power(chi, -1),
            "The Rankin determinant character must be chi^-1.")

    # Independent coordinates eta(m), eta(r), psi(r). In the divisor
    # term of E_2(1,chi), twisting by eta gives
    # eta(m) chi(r) = eta(m/r) psi(r).
    eta_m = (1, 0, 0)
    eta_r = (0, 1, 0)
    psi_r = (0, 0, 1)
    twisted_divisor_term = product(eta_m, power(eta_r, -1), psi_r)
    inducing_pair_divisor_term = product(
        product(eta_m, power(eta_r, -1)), psi_r)
    require(twisted_divisor_term == inducing_pair_divisor_term,
            "The twist does not give the ordered pair (eta,psi).")
    trivial_character = (0, 0, 0)
    twisted_pair = tuple(product(character, eta)
                         for character in (trivial_character, chi))
    reversed_twisted_pair = tuple(product(character, eta)
                                  for character in (chi, trivial_character))
    require(twisted_pair == (eta, psi),
            "The ordered endpoint must twist to (eta,psi).")
    require(reversed_twisted_pair == (psi, eta),
            "A reversed input pair should twist to (psi,eta).")
    require(twisted_pair != reversed_twisted_pair,
            "Generic inducing pairs must retain their ordering.")

    return {
        "character_coordinates": ["eta", "psi", "omega_h"],
        "chi": list(chi),
        "omega_f": list(omega_f),
        "omega_f_eta": list(omega_f_eta),
        "rankin_determinant_character": list(rankin_character),
        "kernel_character": "chi^(-1)",
        "projected_endpoint_central_character": "chi",
        "ordered_endpoint_input": "E_2(1,chi)",
        "twisted_ordered_endpoint": "E_2(eta,psi)",
        "reversed_ordered_endpoint": "E_2(chi,1) would twist to E_2(psi,eta), generally a different form.",
        "divisor_term_identity": "eta(m) chi(r) = eta(m/r) psi(r)",
        "normalization": "Q_chi(w)/u_chi, where u_chi is its actual nonzero unit-ideal coefficient at w=0; no explicit value is asserted.",
        "method": "Exact identities in a free character-exponent lattice; analytic and automorphic identifications are external inputs.",
    }


def verify():
    one = {ZERO: 1}
    A, B, C, D, X = [variable(i) for i in range(5)]
    roots = [multiply(A, C), multiply(A, D),
             multiply(B, C), multiply(B, D)]
    factors = [add(one, scale(multiply(root, X), -1)) for root in roots]

    # h_j(A,B)=(A^(j+1)-B^(j+1))/(A-B), and analogously for C,D.
    # Summing the four geometric series gives the numerator on the left.
    binet_cross_numerator = {}
    for i, sign in enumerate((1, -1, -1, 1)):
        term = roots[i]
        for j, factor in enumerate(factors):
            if j != i:
                term = multiply(term, factor)
        binet_cross_numerator = add(binet_cross_numerator, scale(term, sign))

    abcd = multiply(multiply(A, B), multiply(C, D))
    rankin_numerator = add(one, scale(multiply(abcd, multiply(X, X)), -1))
    binet_denominator = multiply(add(A, scale(B, -1)), add(C, scale(D, -1)))
    claimed_cross_numerator = multiply(binet_denominator, rankin_numerator)
    difference = add(binet_cross_numerator, scale(claimed_cross_numerator, -1))
    require(not difference, "The generic Rankin polynomial identity failed.")

    # The arithmetic Satake products for weights ell+2 and ell are
    # chi^-1 q^(2 ell). X=q^(-(ell+1+w)); therefore ABCD X^2 has
    # q-exponent 2 ell - 2(ell+1+w) = -2-2w.
    abcd_q_exponent = {"ell": 2, "w": 0, "constant": 0}
    x_q_exponent = {"ell": -1, "w": -1, "constant": -1}
    substituted_exponent = {
        key: abcd_q_exponent[key] + 2 * x_q_exponent[key]
        for key in abcd_q_exponent
    }
    require(substituted_exponent == {"ell": 0, "w": -2, "constant": -2},
            "The weight/exponent substitution failed.")
    character_directions = verify_character_directions()

    return {
        "status": "PASS",
        "verification_scope": "Exact local polynomial algebra and formal character-direction identities; no analytic theorem is certified.",
        "python_version": sys.version.split()[0],
        "python_optimization": sys.flags.optimize,
        "coefficient_ring": "Z[A,B,C,D,X]",
        "variables": list(VARIABLES),
        "binet_formula": "h_j(A,B)=(A^(j+1)-B^(j+1))/(A-B)",
        "generic_identity": "sum_j>=0 h_j(A,B) h_j(C,D) X^j = (1-ABCD X^2)/((1-ACX)(1-ADX)(1-BCX)(1-BDX))",
        "method": "Exact polynomial cross multiplication of all four geometric-series denominators and (A-B)(C-D).",
        "cross_numerator": serialize(binet_cross_numerator),
        "claimed_cross_numerator": serialize(claimed_cross_numerator),
        "difference": serialize(difference),
        "repeated_roots": "The cross-multiplied polynomial identity is valid identically; repeated-root generating functions follow by polynomial specialization.",
        "substitution": {
            "ABCD": "chi(p)^(-1) q^(2 ell)",
            "X": "q^(-(ell+1+w))",
            "q_exponent_computation": "2 ell - 2(ell+1+w) = -2-2w",
            "q_exponent_coefficients": substituted_exponent,
            "local_numerator": "1-chi(p)^(-1) q^(-2-2w)",
            "global_correction_character": "chi^(-1)",
            "global_correction_argument": "2+2w",
        },
        "character_directions": character_directions,
        "not_verified": [
            "The identification of Satake products with central characters and weights.",
            "The identification of Mok's projected kernel with the ordered Eisenstein family E_2(1,chi).",
            "The ordinary class-group projection, regularity, and nonzero actual endpoint coefficient u_chi.",
            "Any explicit global scalar, Gauss sum, or different-ideal normalization factor.",
            "Analytic continuation and specialization of the Rankin integral.",
            "Boundary nonvanishing of the Rankin-Selberg L-function.",
        ],
    }


def main():
    if not __debug__ or sys.flags.optimize:
        print("ERROR: optimized Python is not permitted for this verifier.", file=sys.stderr)
        return 2
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the JSON certificate to this path.")
    args = parser.parse_args()
    try:
        result = verify()
    except ArithmeticError as error:
        print("FAIL: " + str(error), file=sys.stderr)
        return 1
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
        print("PASS: exact Rankin local factor and character directions; JSON written to " + str(args.output))
    else:
        print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
