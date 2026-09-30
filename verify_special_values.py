#!/usr/bin/env python3
"""Recompute four dyadic reciprocal constants by exact Dirichlet characters.

Uses generalized Bernoulli sums in Q(zeta_3) and Q(i), without Sage.
The cyclic-cubic identification uses the conductor-discriminant theorem and
uniqueness of the order-three characters modulo the prime conductor.
The quartic field is identified directly as Q(zeta_15)^+ by a polynomial identity.
Artin factorization and the generalized Bernoulli special-value formula are
external mathematical theorems, explicitly not proved by this arithmetic check.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit('Do not use optimized Python; assertions are required.')

script_directory = str(Path(__file__).resolve().parent)
if script_directory not in sys.path:
    sys.path.insert(0, script_directory)
from verify_dyadic_bounds import Field, FIELDS


def primitive_logs(p, g):
    result = {pow(g, j, p): j for j in range(p-1)}
    assert len(result) == p-1
    return result


def bernoulli_pair(conductor, values, order):
    """Coordinates in basis 1,zeta_3 or 1,i."""
    sums = [[0, 0] for _ in range(3)]
    for a, j in values.items():
        if order == 3:
            c = [(1, 0), (0, 1), (-1, -1)][j % 3]
        else:
            c = [(1, 0), (0, 1), (-1, 0), (0, -1)][j % 4]
        for k in range(3):
            for i in range(2):
                sums[k][i] += c[i]*a**k
    assert sums[0] == sums[1] == [0, 0]
    return [-Q(c, 2*conductor) for c in sums[2]]


def Laurent_power(j):
    result = {0: 1}
    for _ in range(j):
        new = {}
        for e, c in result.items():
            for shift in (-1, 1):
                new[e+shift] = new.get(e+shift, 0)+c
        result = new
    return result


def verify():
    rows = []
    for d, p, generator, expected in [(49, 7, 3, -168), (169, 13, 2, -24), (361, 19, 2, -8)]:
        assert d == p*p
        Field(FIELDS[d][0], d)  # irreducibility and square discriminant
        logs = primitive_logs(p, generator)
        char = {a: j % 3 for a, j in logs.items()}
        assert char[p-1] == 0
        for a in char:
            for b in char:
                assert char[a*b % p] == (char[a]+char[b]) % 3
        L = bernoulli_pair(p, char, 3)
        norm = L[0]**2-L[0]*L[1]+L[1]**2
        zeta = -norm/12
        alpha = 8/zeta
        assert alpha == expected
        rows.append({'discriminant': d, 'primitive_conductor': p,
                     'character_exponents': char, 'L_minus_one_coordinates': list(map(str, L)),
                     'zeta_minus_one': str(zeta), 'alpha': str(alpha)})
    f = FIELDS[1125][0]
    Field(f, 1125)
    # z^4 f(z+z^-1)=Phi_15(z), checked coefficient by coefficient.
    poly = {}
    for j, c in enumerate(f):
        for e, v in Laurent_power(j).items():
            poly[e+4] = poly.get(e+4, 0)+c*v
    phi15 = [1, -1, 0, 1, -1, 1, 0, -1, 1]
    assert [poly.get(j, 0) for j in range(9)] == phi15
    assert not any(v for k, v in poly.items() if k < 0 or k > 8)
    char = {((-1)**sgn*pow(2, j, 15)) % 15: j for sgn in (0, 1) for j in range(4)}
    assert len(char) == 8 and char[14] == 0
    for a in char:
        for b in char:
            assert char[a*b % 15] == (char[a]+char[b]) % 4
    # The square character is induced from the primitive quadratic character mod5.
    squares = {a*a % 5 for a in range(1, 5)}
    chi5 = {a: 1 if a in squares else -1 for a in range(1, 5)}
    for a, j in char.items():
        assert (-1)**j == chi5[a % 5]
    L5 = -Q(sum(chi5[a]*a*a for a in chi5), 10)
    L = bernoulli_pair(15, char, 4)
    zeta = -Q(1, 12)*L5*(L[0]**2+L[1]**2)
    alpha = 16/zeta
    assert zeta == Q(4, 15) and alpha == 60
    rows.append({'discriminant': 1125, 'cyclotomic_conductor': 15,
                 'primitive_conductors': [1, 5, 15, 15],
                 'real_cyclotomic_polynomial_identity': phi15,
                 'character_exponents': char, 'L_minus_one_coordinates': list(map(str, L)),
                 'L_minus_one_quadratic_primitive': str(L5),
                 'zeta_minus_one': str(zeta), 'alpha': str(alpha)})
    return {'status': 'passed', 'arithmetic': 'integer/rational, Python standard library',
            'external_dependencies': ['Artin factorization', 'conductor-discriminant theorem',
                                      'generalized Bernoulli formula'], 'fields': rows}


def main():
    if not __debug__:
        raise SystemExit('Assertions must be enabled: do not use -O/-OO.')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    result = verify()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: exact reciprocal constants -168,-24,-8,60')


if __name__ == '__main__':
    main()
