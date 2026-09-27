#!/usr/bin/env python3
"""Free, rigorous S_2=0 certificate for the quartic field of discriminant725.

Eichler's class-number formula: Kirschmer--Voight, Algorithmic enumeration
of ideal classes for quaternion orders, Proposition5.1 and Lemma5.3.
The mass is rigorously bounded, not numerically rationally reconstructed.
The unique integer in the resulting class-number interval is1.
"""
from sage.all import *
from sage.version import version as sage_version
import argparse
import json

proof.all(True)


def run():
    x = polygen(QQ)
    f = x**4-x**3-3*x*x+x+1
    F = NumberField(f, 'a')
    assert F.discriminant() == 725
    assert F.class_number(proof=True) == F.narrow_class_group().order() == 1
    t = F['t'].gen()
    # For K/F quadratic, phi(2q)|8; phi(m)^2>=m/2 implies2q<=128.
    candidates = [q for q in range(2,65) if 8 % euler_phi(2*q) == 0]
    cm_rows = []
    correction = QQ(0)
    excluded = []
    for q in candidates:
        cyclo = F['t'](cyclotomic_polynomial(2*q))
        factors = list(cyclo.factor())
        quadratics = [g for g,e in factors if g.degree() == 2]
        if not quadratics:
            excluded.append({'q': q, 'factor_degrees': [int(g.degree()) for g,e in factors]})
            continue
        assert q in (2,3,5)
        g = quadratics[0]
        E = F.extension(g, 'b')
        # All conjugates of this cyclotomic generator lie in E. Thus all
        # degree-two factors describe the same quadratic extension of F.
        assert len(E['s'](cyclo).roots()) == euler_phi(2*q)
        A = E.absolute_field('c')
        rel_discriminant = E.relative_discriminant()
        polynomial_discriminant = F.ideal(g.discriminant())
        assert rel_discriminant == polynomial_discriminant
        # Discriminant equality means O_F[zeta_(2q)] is maximal.
        cm_h = A.class_number(proof=True)
        assert cm_h == 1
        roots_of_unity = A.number_of_roots_of_unity()
        assert roots_of_unity == 2*q
        correction += QQ(cm_h)/2*(1-QQ(1)/q)
        cm_rows.append({
            'q': int(q), 'relative_polynomial': str(g),
            'absolute_polynomial': str(A.defining_polynomial()),
            'absolute_discriminant': int(A.discriminant()),
            'relative_discriminant_norm': int(rel_discriminant.norm()),
            'polynomial_discriminant_norm': int(polynomial_discriminant.norm()),
            'cyclotomic_order_is_maximal': True,
            'all_cyclotomic_conjugates_in_same_extension': True,
            'number_of_roots_of_unity': int(roots_of_unity),
            'class_number': int(cm_h), 'unit_index_mod_base_units': int(q),
        })
    assert [r['q'] for r in cm_rows] == [2,3,5]
    assert correction == QQ(59)/60
    # A crude Euler bound already suffices. zeta_F(2) <= zeta(2)^4,
    # and zeta(2) < 2 by sum_{m>=2} m^-2 < integral_1^infty t^-2 dt.
    # Functional equation yields mass=2 D^(3/2)zeta_F(2)/(4*pi^2)^4.
    # pi>3 gives the fully rational/algebraic strict upper bound below.
    mass_upper = AA(2*725)*AA(725).sqrt()*16/36**4
    assert 0 < mass_upper < 1
    class_number_lower = correction
    class_number_upper = AA(correction)+mass_upper
    assert 0 < class_number_lower < 1 < class_number_upper < 2
    # Strict interval (59/60,59/60+mass_upper) has the unique integer1.
    return {
        'sage_version': sage_version, 'proof_all': True,
        'field_polynomial': str(f), 'field_discriminant': 725,
        'ordinary_class_number': 1, 'narrow_class_number': 1,
        'elliptic_order_candidates': candidates, 'excluded_orders': excluded,
        'CM_certificates': cm_rows, 'elliptic_correction': str(correction),
        'mass_lower_strict': '0',
        'mass_upper_expression': '2*725*sqrt(725)*16/36^4',
        'mass_upper_interval': str(RealIntervalField(128)(mass_upper)),
        'class_number_lower_strict': str(class_number_lower),
        'class_number_upper_interval': str(RealIntervalField(128)(class_number_upper)),
        'unique_integer_class_number': 1,
        'mass_deduced_from_class_number_formula': '1/60',
        'zeta_minus_one_deduced_from_mass': '2/15',
        'full_level_parallel_weight_two_cusp_dimension': 0,
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output')
    args = p.parse_args()
    result = json.dumps(run(), indent=2)
    if args.output:
        with open(args.output, 'w') as handle:
            handle.write(result+'\n')
    else:
        print(result)
