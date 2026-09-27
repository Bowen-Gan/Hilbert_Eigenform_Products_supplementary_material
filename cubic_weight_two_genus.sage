#!/usr/bin/env python3
"""Free, exact weight-two vanishing certificate for D=169 and D=361.

Uses Shimizu's area formula and the elliptic-cycle formula in J. Voight,
Shimura curves of genus at most two, equations (1),(2), Lemma 2.1 and
Proposition 2.3(a).  The mathematical reduction is in the companion note.
No Hilbert-modular-form implementation, database, or Magma is used.

Run: sage -python cubic_weight_two_genus.sage --output cubic_weight_two_genus.json
"""
from sage.all import *
from sage.version import version as sage_version
import argparse
import json

proof.all(True)
x = polygen(QQ)


def exact_cubic_zeta(K, conductor):
    C = CyclotomicField(3)
    Z = CyclotomicField(conductor)
    z = Z.gen()
    for chi in DirichletGroup(conductor, C):
        if chi.order() != 3 or chi(-1) != 1:
            continue
        period = sum(z**a for a in range(1, conductor) if chi(a) == 1)
        f = period.minpoly()
        if f.degree() != 3 or not NumberField(f, 'b').is_isomorphic(K):
            continue
        L_values = []
        value = C(1)
        for j in range(3):
            c = (chi**j).primitive_character()
            m = c.conductor()
            L = -sum(c(a)*(QQ(a*a)/m-a+QQ(m)/6)
                     for a in range(1, m+1))/2
            assert L == -c.bernoulli(2)/2
            value *= L
            L_values.append(str(L))
        return QQ(value), {
            'conductor': conductor, 'period_minpoly': str(f),
            'field_isomorphism_verified': True, 'primitive_L_values': L_values,
            'coefficient_field': str(C.defining_polynomial()),
        }
    raise AssertionError('No exact Gaussian-period field identification')


def run():
    rows = []
    for D, f, conductor in [
        (169, x**3-x**2-4*x-1, 13),
        (361, x**3-x**2-6*x+7, 19),
    ]:
        F = NumberField(f, 'a')
        assert F.discriminant() == D
        h = F.class_number(proof=True)
        hp = F.narrow_class_group().order()
        assert h == hp == 1
        zeta, zeta_certificate = exact_cubic_zeta(F, conductor)
        assert zeta == {169: QQ(-1)/3, 361: QQ(-1)}[D]
        # If F(zeta_(2q))/F is quadratic then phi(2q) divides 6.
        # phi(m)^2 >= m/2 implies m<=72.  Thus this finite loop is complete.
        possible_q = [q for q in range(2,37) if 6 % euler_phi(2*q) == 0]
        assert possible_q == [2,3,7,9]
        excluded_real_subfields = []
        for q in (7,9):
            Z = CyclotomicField(2*q)
            eta = Z.gen()+Z.gen()**(-1)
            Kplus = NumberField(eta.minpoly(), 'r')
            assert Kplus.degree() == 3
            assert Kplus.discriminant() != D
            excluded_real_subfields.append({
                'q': q, 'real_cyclotomic_discriminant': int(Kplus.discriminant())})
        t = F['t'].gen()
        cm_rows = []
        counts = {}
        for q, cm_poly, rational_discriminant in [(2,t*t+1,-4),(3,t*t+t+1,-3)]:
            E = F.extension(cm_poly, 'u').absolute_field('c')
            assert E.degree() == 6
            # Equality with the tensor-order discriminant proves maximality
            # of O_F[zeta_(2q)]. Therefore it is the only relevant order.
            tensor_discriminant = ZZ(D)**2 * ZZ(rational_discriminant)**3
            assert E.discriminant() == tensor_discriminant
            cm_h = E.class_number(proof=True)
            assert cm_h == {(169,2):3,(169,3):1,(361,2):1,(361,3):3}[(D,q)]
            counts[q] = cm_h
            cm_rows.append({
                'elliptic_order': q,
                'absolute_polynomial': str(E.defining_polynomial()),
                'absolute_discriminant': int(E.discriminant()),
                'tensor_order_discriminant': int(tensor_discriminant),
                'cyclotomic_order_is_maximal': True,
                'class_number': int(cm_h), 'norm_unit_index': 1,
                'elliptic_cycles': int(cm_h),
            })
        area = abs(zeta)/2
        genus = 1+abs(zeta)/4-QQ(counts[2])/4-QQ(counts[3])/3
        assert genus == 0
        assert area == 2*genus-2+counts[2]*QQ(1)/2+counts[3]*QQ(2)/3
        rows.append({
            'discriminant': D, 'polynomial': str(f),
            'class_number': int(h), 'narrow_class_number': int(hp),
            'zeta_minus_one': str(zeta), 'zeta_certificate': zeta_certificate,
            'possible_elliptic_orders': possible_q,
            'excluded_cubic_real_subfields': excluded_real_subfields,
            'CM_certificates': cm_rows, 'normalized_area': str(area),
            'genus': int(genus), 'full_level_parallel_weight_two_cusp_dimension': 0,
        })
    return {'sage_version': sage_version, 'proof_all': True,
            'method': 'Shimizu area, exact CM class numbers, Riemann-Hurwitz, Jacquet-Langlands',
            'rows': rows}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output')
    args = p.parse_args()
    output = json.dumps(run(), indent=2)
    if args.output:
        with open(args.output, 'w') as handle:
            handle.write(output+'\n')
    else:
        print(output)
