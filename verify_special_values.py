#!/usr/bin/env python3
"""Recompute the manuscript's finite Dirichlet special-value tables exactly.

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
from math import factorial
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit('Do not use optimized Python; assertions are required.')

script_directory = str(Path(__file__).resolve().parent)
if script_directory not in sys.path:
    sys.path.insert(0, script_directory)
from verify_dyadic_bounds import Field, FIELDS
from rq_certificate import kronecker


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


def critical_value_pair(conductor, values, order, weight):
    """L(1-weight,chi) in Q(zeta_3) or Q(i), with full Bernoulli sums.

    The coefficient lists below are B_3(X) and B_4(X), respectively.
    No vanishing first moment is assumed for an odd character.
    """
    bernoulli = {3: [Q(0), Q(1, 2), Q(-3, 2), Q(1)],
                 4: [Q(-1, 30), Q(0), Q(1), Q(-2), Q(1)]}[weight]
    result = [Q(0), Q(0)]
    roots = ([(1, 0), (0, 1), (-1, -1)] if order == 3
             else [(1, 0), (0, 1), (-1, 0), (0, -1)])
    assert order in (3, 4)
    for a, j in values.items():
        b = sum(c*Q(a, conductor)**k for k, c in enumerate(bernoulli))
        multiplier = -Q(conductor**(weight-1), weight)*b
        for i, c in enumerate(roots[j % order]):
            result[i] += c*multiplier
    return result


def quadratic_critical_value(discriminant, weight):
    conductor = abs(discriminant)
    bernoulli = {3: [Q(0), Q(1, 2), Q(-3, 2), Q(1)],
                 4: [Q(-1, 30), Q(0), Q(1), Q(-2), Q(1)]}[weight]
    total = sum(kronecker(discriminant, a)*
                sum(c*Q(a, conductor)**k for k, c in enumerate(bernoulli))
                for a in range(1, conductor+1))
    return -Q(conductor**(weight-1), weight)*total


def verify_rational_cutoffs():
    """Integer/rational comparisons after the stated pi and zeta bounds."""
    pi_upper = Q(22, 7)
    comparisons = []
    for D, weight, bound in [(28, 4, 1), (21, 5, 1), (12, 4, 14)]:
        V_squared = (Q(4**(2*weight)*D,
                         8100**2*D**(2*weight)*factorial(weight-1)**4)*
                     pi_upper**(4*weight+16))
        assert V_squared < bound**2
        comparisons.append({'comparison': f'V({D},{weight})<{bound}',
                            'squared_upper': str(V_squared)})
    # b_k(delta)^2 avoids fractional powers entirely.
    for weight, delta, zeta_denominator, bound in [
        (4, Q(51, 10), 90, Q(19, 20)),
        (6, Q(18, 5), 945, Q(1, 2)),
    ]:
        zeta_upper = pi_upper**weight/zeta_denominator
        b_squared = (2*pi_upper*zeta_upper**2/factorial(weight-1)**2*
                     (2*pi_upper/delta)**(2*weight-1))
        assert b_squared < bound**2
        assert 2*pi_upper/(delta*weight) < 1
        comparisons.append({'comparison': f'b_{weight}({delta})<{bound}',
                            'squared_upper': str(b_squared)})
    zeta3_upper = sum(Q(1, j**3) for j in range(1, 26))+Q(1, 1250)
    assert zeta3_upper < Q(12021, 10000)
    cutoff_base = 4*Q(355, 113)**3*zeta3_upper
    assert cutoff_base**4 < 55**5
    assert cutoff_base**8 < 3003**5
    assert cutoff_base**2 < Q(7403, 1000)**5
    for D, bound in [(28, 6), (24, 9)]:
        squared_upper = Q(25**2)*pi_upper**12/D**5
        assert squared_upper < bound**2
        comparisons.append({'comparison': f'25*pi^6/{D}^(5/2)<{bound}',
                            'squared_upper': str(squared_upper)})
    return {'external_pi_bounds': ['pi<22/7', 'pi<355/113'],
            'zeta_three_upper_by_integral_test': str(zeta3_upper),
            'delta_three_checks': ['delta_3^2<55', 'delta_3^4<3003', 'delta_3<7.403'],
            'exact_squared_comparisons': comparisons}


def verify_higher_eisenstein_weights():
    """Independent Bernoulli checks for Eisenstein weights three and four."""
    cubic_rows = []
    for discriminant, conductor, group_order, expected_product, expected_zeta in [
        (49, 7, 6, Q(316, 7), Q(79, 210)),
        (81, 9, 6, Q(796, 3), Q(199, 90)),
    ]:
        polynomial = FIELDS[49][0] if discriminant == 49 else [-1, -3, 0, 1]
        Field(polynomial, discriminant)
        generator = 3 if conductor == 7 else 2
        logs = {pow(generator, j, conductor): j for j in range(group_order)}
        assert len(logs) == group_order
        character = {a: j % 3 for a, j in logs.items()}
        assert character[conductor-1] == 0
        if conductor == 9:
            # A character induced modulo 3 would take the same value at
            # 1 and 4; the cubic character does not.
            assert character[1] != character[4]
        for a in character:
            for b in character:
                assert character[a*b % conductor] == (character[a]+character[b]) % 3
        L = critical_value_pair(conductor, character, 3, 4)
        norm = L[0]**2-L[0]*L[1]+L[1]**2
        zeta = norm/120
        assert norm == expected_product and zeta == expected_zeta
        alpha = 8/zeta
        assert alpha.denominator != 1
        cubic_rows.append({'discriminant': discriminant, 'conductor': conductor,
                           'L_minus_three_coordinates': list(map(str, L)),
                           'conjugate_product': str(norm), 'zeta_minus_three': str(zeta),
                           'alpha': str(alpha)})
    L21 = quadratic_critical_value(21, 4)
    assert L21 == 308 and -4*L21 == -1232

    quadratic_rows = []
    for D, d1, d2, expected, expected_alpha in [
        (12, -3, -4, Q(1, 9), Q(36)),
        (21, -3, -7, Q(32, 63), Q(63, 8)),
        (24, -3, -8, Q(2, 3), Q(6)),
        (28, -4, -7, Q(8, 7), Q(7, 2)),
        (33, -3, -11, Q(4, 3), Q(3)),
        (44, -4, -11, Q(3), Q(4, 3)),
    ]:
        factors = [quadratic_critical_value(d, 3) for d in (d1, d2)]
        value = factors[0]*factors[1]
        assert d1*d2 == D and value == expected and 4/value == expected_alpha
        quadratic_rows.append({'discriminant': D, 'genus_discriminants': [d1, d2],
                               'primitive_factors': list(map(str, factors)),
                               'L_minus_two': str(value), 'alpha': str(4/value)})

    # The odd quartic characters modulo 5 are primitive and conjugate.
    character5 = {1: 0, 2: 1, 4: 2, 3: 3}
    L5 = critical_value_pair(5, character5, 4, 3)
    product5 = L5[0]**2+L5[1]**2
    assert product5 == Q(4, 5)
    quartic_rows = []
    for conductor, D, negative_discriminants, extra_factor, expected in [
        (15, 1125, [-3, -15], product5, Q(128, 45)),
        (20, 2000, [-4, -20], product5, Q(12)),
        (24, 2304, [-4, -8, -3, -24], Q(1), Q(46, 3)),
    ]:
        factors = [quadratic_critical_value(d, 3) for d in negative_discriminants]
        value = extra_factor
        for factor in factors:
            value *= factor
        assert value == expected and (16/value).denominator != 1
        quartic_rows.append({'discriminant': D, 'cyclotomic_conductor': conductor,
                             'quadratic_primitive_discriminants': negative_discriminants,
                             'quadratic_primitive_factors': list(map(str, factors)),
                             'conjugate_conductor_five_product': str(extra_factor),
                             'L_minus_two': str(value), 'alpha': str(16/value)})
    return {'weight_four_cubic_fields': cubic_rows,
            'weight_four_D21': {'B_4': '-1232', 'L_minus_three': str(L21),
                                'zeta_minus_three': str(L21/120), 'alpha': str(4/(L21/120))},
            'weight_three_quadratic_fields': quadratic_rows,
            'weight_three_quartic_fields': quartic_rows,
            'rational_numerical_cutoffs': verify_rational_cutoffs()}


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
                                      'generalized Bernoulli formula'], 'fields': rows,
            'higher_eisenstein_weight_checks': verify_higher_eisenstein_weights()}


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
    print('PASS: exact dyadic reciprocal constants and Eisenstein-weight-three/four tables')


if __name__ == '__main__':
    main()
