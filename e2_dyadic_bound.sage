#!/usr/bin/env python3
"""Case 11: exact dyadic certificates for the surviving small inert fields.

Run: sage -python e2_dyadic_bound.sage --output e2_dyadic_output.json
or use a Python installation providing sage.all.

No special value is recovered from numerical approximation.  The five
abelian special values are generalized Bernoulli sums after an exact Gaussian
period field identification.  For D=725,5125 a finite Euler product supplies a
rigorous upper bound for |2^4/zeta_F(-1)| via the functional equation.
"""
from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
from sage.version import version as sage_version
import argparse
import json
import sys

proof.all(True)
x = polygen(QQ)
ROWS = [
    (49, x**3-x**2-2*x+1, 7),
    (169, x**3-x**2-4*x-1, 13),
    (361, x**3-x**2-6*x+7, 19),
    (725, x**4-x**3-3*x**2+x+1, None),
    (1125, x**4-x**3-4*x**2+4*x+1, 15),
    (5125, x**4-2*x**3-6*x**2+7*x+11, None),
    (6125, x**4-x**3-9*x**2+9*x+11, 35),
]
RIF256 = RealIntervalField(256)


def abelian_special_value(K, conductor):
    """Identify K as the fixed field of ker(chi), then use all its characters."""
    n = int(K.degree())
    C = CyclotomicField(euler_phi(conductor))
    Z = CyclotomicField(conductor)
    z = Z.gen()
    for chi in DirichletGroup(conductor, C):
        if chi.order() != n or chi(-1) != 1 or chi.conductor() != conductor:
            continue
        kernel = [a for a in range(1, conductor) if chi(a) == 1]
        period_polynomial = sum(z**a for a in kernel).minpoly()
        if period_polynomial.degree() != n:
            continue
        period_field = NumberField(period_polynomial, 'b')
        if not period_field.is_isomorphic(K):
            continue
        character_rows = []
        zeta_value = C(1)
        for j in range(n):
            primitive = (chi**j).primitive_character()
            f = primitive.conductor()
            # B_{2,chi}=f sum_{a=1}^f chi(a) B_2(a/f).
            bernoulli_sum = sum(
                primitive(a)*(QQ(a*a)/f-a+QQ(f)/6)
                for a in range(1, f+1)
            )
            assert bernoulli_sum == primitive.bernoulli(2)
            L = -bernoulli_sum/2
            zeta_value *= L
            character_rows.append({
                'power': j, 'primitive_conductor': int(f),
                'values_on_residues': [str(primitive(a)) for a in range(1, f+1)],
                'L_minus_one': str(L),
            })
        value = QQ(zeta_value)
        return value, {
            'ambient_conductor': conductor,
            'coefficient_field': str(C.defining_polynomial()),
            'kernel_residues': kernel,
            'period_minpoly': str(period_polynomial),
            'field_isomorphism_checked': True,
            'characters': character_rows,
            'zeta_minus_one': str(value),
            'alpha': str(QQ(2**n)/value),
        }
    raise AssertionError('No exact character-field identification found')


def euler_alpha_upper(K, cutoff=199):
    """Functional equation: |alpha|=(4*pi^2)^n/(D^(3/2)*zeta_F(2))."""
    n = int(K.degree())
    D = ZZ(K.discriminant())
    product = QQ(1)
    factors = []
    for p in prime_range(cutoff+1):
        local = QQ(1)
        norms = []
        for P, e in K.ideal(p).factor():
            N = ZZ(P.norm())
            local /= 1-QQ(1)/N**2
            norms.append({'norm': int(N), 'ramification_index': int(e)})
        product *= local
        factors.append({'p': int(p), 'primes': norms, 'factor': str(local)})
    # MPFI gives an outward-rounded upper endpoint for pi. All following
    # operations and the final rational ceiling are exact real algebraic ones.
    pi_upper = QQ(RIF256.pi().upper())
    exact_upper = AA((4*pi_upper**2)**n)/(AA(D)*AA(D).sqrt()*product)
    rational_upper = QQ(ceil(exact_upper*10**9))/10**9
    assert AA(rational_upper) >= exact_upper
    return rational_upper, {
        'cutoff': cutoff, 'Euler_factors': factors,
        'Euler_product': str(product), 'pi_upper': str(pi_upper),
        'alpha_absolute_upper': str(rational_upper),
        'method': 'finite Euler lower product and functional equation',
    }


def divisor_data(I):
    factors = list(I.factor())
    tau = prod(e+1 for P, e in factors)
    sigma = prod(sum(P.norm()**j for j in range(e+1)) for P, e in factors)
    return ZZ(tau), ZZ(sigma), [(int(P.norm()), int(e)) for P, e in factors]


def positive_decompositions(K):
    basis = K.integral_basis()
    n = int(K.degree())
    gram = matrix(ZZ, [[(u*v).trace() for v in basis] for u in basis])
    lattice = IntegralLattice(gram)
    # If 0 << u << 4 then Tr(u^2)<16n. short_vectors is exhaustive
    # below this strict bound and includes both signs.
    bound = 16*n
    terms = []
    examined = 0
    for layer in lattice.short_vectors(bound, up_to_sign_flag=False):
        for coordinates in layer:
            examined += 1
            u = sum((ZZ(c)*b for c, b in zip(coordinates, basis)), K(0))
            v = K(4)-u
            if not u.is_totally_positive() or not v.is_totally_positive():
                continue
            tx, sx, fx = divisor_data(K.ideal(u))
            ty, sy, fy = divisor_data(K.ideal(v))
            Nv = ZZ(v.norm())
            assert 0 < Nv < 4**n
            terms.append({
                'coordinates': [int(c) for c in coordinates],
                'x': str(u), 'y': str(v),
                'norm_x': int(u.norm()), 'norm_y': int(Nv),
                'sigma_one_x': int(sx), 'divisor_count_y': int(ty),
                'ideal_x_factor_norms': fx, 'ideal_y_factor_norms': fy,
            })
    terms.sort(key=lambda t: t['coordinates'])
    return terms, {
        'integral_basis': [str(b) for b in basis],
        'trace_gram_matrix': [[int(c) for c in row] for row in gram.rows()],
        'strict_squared_length_bound': bound,
        'number_of_short_vectors_including_zero': examined,
        'number_of_positive_decompositions': len(terms),
    }


def normalized_right(terms, q, alpha_bound, t):
    convolution = sum(AA(r['sigma_one_x']*r['divisor_count_y']) *
                      (AA(r['norm_y']).sqrt()/q)**t for r in terms)
    return AA(alpha_bound)*(convolution+4*AA(q).sqrt()**(-t)
                             +AA(alpha_bound)*AA(q)**(-t))


def run():
    output = {'sage_version': sage_version, 'proof_all': True, 'fields': []}
    expected_alpha = {49: -168, 169: -24, 361: -8,
                      1125: 60, 6125: QQ(60)/13}
    for D, f, conductor in ROWS:
        K = NumberField(f, 'a')
        n = int(K.degree())
        assert K.discriminant() == D
        assert K.class_number(proof=True) == 1
        narrow_class_number = int(K.narrow_class_group().order())
        if D in (49,169,361,725):
            assert narrow_class_number == 1
        factor_two = list(K.ideal(2).factor())
        assert len(factor_two) == 1
        assert factor_two[0][1] == 1 and factor_two[0][0].norm() == 2**n
        row = {'discriminant': D, 'polynomial': str(f), 'degree': n,
               'ordinary_class_number': 1, 'narrow_class_number': narrow_class_number,
               'dyadic_norm': int(2**n)}
        if conductor:
            zeta, value_certificate = abelian_special_value(K, conductor)
            alpha = QQ(2**n)/zeta
            assert alpha == expected_alpha[D]
            row['special_value_certificate'] = value_certificate
            alpha_bound = abs(alpha)
            if alpha.denominator() != 1:
                row['exclusion'] = 'nonintegral reciprocal constant term'
                output['fields'].append(row)
                continue
        else:
            alpha_bound, value_certificate = euler_alpha_upper(K)
            row['alpha_upper_certificate'] = value_certificate
        terms, lattice_certificate = positive_decompositions(K)
        row['lattice_certificate'] = lattice_certificate
        row['decompositions'] = terms
        q = ZZ(2**n)
        left = AA(q*q-1)
        comparisons = []
        threshold = None
        for t in range(1, 101):
            rhs = normalized_right(terms, q, alpha_bound, t)
            smaller = bool(rhs < left)
            comparisons.append({'t': t, 'source_weight': t+1,
                                'left': str(left),
                                'right_interval': str(RIF256(rhs)),
                                'right_strictly_less': smaller})
            if smaller:
                threshold = t+1
                break
        assert threshold is not None
        assert threshold == {49: 7, 169: 4, 361: 3, 725: 4, 1125: 4, 5125: 2}[D]
        row['alpha_absolute_bound_used'] = str(alpha_bound)
        row['comparisons'] = comparisons
        row['all_source_weights_at_least_excluded'] = threshold
        # Every norm_y<q^2 and q>1: every term in the normalized RHS
        # strictly decreases with t. This proves all subsequent weights.
        row['monotonicity_bases_less_than_one_checked'] = True
        output['fields'].append(row)
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default=None)
    args = parser.parse_args()
    result = run()
    encoded = json.dumps(result, indent=2)
    if args.output:
        with open(args.output, 'w') as handle:
            handle.write(encoded+'\n')
    else:
        print(encoded)
    for row in result['fields']:
        print('D=%s: %s' % (row['discriminant'], row.get('exclusion',
              'source weights >= %s excluded' % row.get('all_source_weights_at_least_excluded'))),
              file=sys.stderr)
