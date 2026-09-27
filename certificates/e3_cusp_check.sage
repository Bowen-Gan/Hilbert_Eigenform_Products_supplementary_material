"""Case 10: exact number-field and primitive Dirichlet L-value certificates.
Run: sage -python e3_cusp_check.sage
This supplements only the finite arithmetic used in the manuscript.
"""

from sage.all import *
from sage.rings.number_field.totallyreal_rel import (
    enumerate_totallyreal_fields_all,
)

proof.all(True)
x = polygen(QQ)


def exact_L_one_minus_k(chi, k):
    """Return L(1-k, chi_primitive) = -B_{k,chi}/k exactly."""
    primitive = chi.primitive_character()
    return -primitive.bernoulli(k) / k


# --------------------------------------------------------------------------
# Analytic cutoff
# --------------------------------------------------------------------------

RBF = RealBallField(200)
R3 = (4 * RBF.pi() ** 3 * RBF.zeta(3)) ** (RBF(2) / 5)
print("R3 = (4*pi^3*zeta(3))^(2/5) =", R3)
print("R3^2 =", R3 ** 2)
print("R3^4 =", R3 ** 4)
assert R3 < RBF(15) / 2
assert R3 ** 2 < 55
assert R3 ** 4 < 3003


# --------------------------------------------------------------------------
# Totally real quartic fields below the forced discriminant cutoff
# --------------------------------------------------------------------------

quartic_fields = enumerate_totallyreal_fields_all(
    4, 3003, return_pari_objects=False
)
expected_quartic_fields = [
    [725, x**4 - x**3 - 3*x**2 + x + 1],
    [1125, x**4 - x**3 - 4*x**2 + 4*x + 1],
    [1600, x**4 - 6*x**2 + 4],
    [1957, x**4 - 4*x**2 - x + 1],
    [2000, x**4 - 5*x**2 + 5],
    [2048, x**4 - 4*x**2 + 2],
    [2225, x**4 - x**3 - 5*x**2 + 2*x + 4],
    [2304, x**4 - 4*x**2 + 1],
    [2525, x**4 - 2*x**3 - 4*x**2 + 5*x + 5],
    [2624, x**4 - 2*x**3 - 3*x**2 + 2*x + 1],
    [2777, x**4 - x**3 - 4*x**2 + x + 2],
]
assert quartic_fields == expected_quartic_fields

print("\nTotally real quartic fields through discriminant 3003")
quartic_without_negative_norm_unit = []
quartic_by_discriminant = {}
for discriminant, polynomial in quartic_fields:
    K = NumberField(polynomial, "a")
    quartic_by_discriminant[ZZ(discriminant)] = K
    unit_norms = [ZZ(u.norm()) for u in K.units()]
    has_negative_norm_unit = -1 in unit_norms
    print(
        discriminant,
        polynomial,
        "unit norms =",
        unit_norms,
        "negative norm =",
        has_negative_norm_unit,
    )
    if not has_negative_norm_unit:
        quartic_without_negative_norm_unit.append(ZZ(discriminant))

assert quartic_without_negative_norm_unit == [1125, 2000, 2304]


# The three survivors are maximal real cyclotomic subfields.  Their full
# cyclotomic extensions have relative finite discriminant one and therefore
# realize the unique nontrivial narrow class character.
quartic_cyclotomic_data = {
    1125: (15, QQ(128) / 45, QQ(45) / 8),
    2000: (20, QQ(12), QQ(4) / 3),
    2304: (24, QQ(46) / 3, QQ(24) / 23),
}

print("\nQuartic class groups and exact odd-character products")
for discriminant, (conductor, expected_L, expected_alpha) in (
    quartic_cyclotomic_data.items()
):
    F = quartic_by_discriminant[ZZ(discriminant)]
    K = CyclotomicField(conductor)
    F_cyclotomic = K.maximal_totally_real_subfield()[0]

    assert F.is_isomorphic(F_cyclotomic)
    assert ZZ(K.discriminant()).abs() == ZZ(discriminant) ** 2
    assert F.class_group().order() == 1
    assert F.narrow_class_group().order() == 2
    assert all(ZZ(u.norm()) == 1 for u in F.units())

    odd_characters = [chi for chi in DirichletGroup(conductor) if chi(-1) == -1]
    odd_values = [exact_L_one_minus_k(chi, 3) for chi in odd_characters]
    print("primitive conductors =", [chi.primitive_character().conductor() for chi in odd_characters])
    L_value = prod(odd_values)
    alpha = QQ(16) / L_value

    print(
        "D =",
        discriminant,
        "m =",
        conductor,
        "odd L(-2) factors =",
        odd_values,
        "L_F(-2,epsilon) =",
        L_value,
        "alpha =",
        alpha,
    )

    assert len(odd_characters) == 4
    assert L_value == expected_L
    assert alpha == expected_alpha
    assert alpha not in ZZ


# --------------------------------------------------------------------------
# Real quadratic fields below the forced cutoff
# --------------------------------------------------------------------------

fundamental_discriminants_below_55 = [
    D for D in range(2, 55) if ZZ(D).is_fundamental_discriminant()
]
expected_fundamental_discriminants = [
    5,
    8,
    12,
    13,
    17,
    21,
    24,
    28,
    29,
    33,
    37,
    40,
    41,
    44,
    53,
]
assert fundamental_discriminants_below_55 == expected_fundamental_discriminants

print("\nReal quadratic fields below discriminant 55")
quadratic_without_negative_norm_unit = []
quadratic_by_discriminant = {}
for discriminant in fundamental_discriminants_below_55:
    if discriminant % 4 == 1:
        polynomial = x**2 - x + (1 - discriminant) // 4
    else:
        polynomial = x**2 - discriminant // 4
    F = NumberField(polynomial, "b")
    quadratic_by_discriminant[ZZ(discriminant)] = F
    assert F.discriminant() == discriminant
    unit_norms = [ZZ(u.norm()) for u in F.units()]
    has_negative_norm_unit = -1 in unit_norms
    print(
        discriminant,
        polynomial,
        "unit norms =",
        unit_norms,
        "negative norm =",
        has_negative_norm_unit,
    )
    if not has_negative_norm_unit:
        quadratic_without_negative_norm_unit.append(ZZ(discriminant))

assert quadratic_without_negative_norm_unit == [12, 21, 24, 28, 33, 44]


# Each quadratic survivor has ordinary class number one and narrow class
# number two.  Its unique nontrivial narrow character is the genus character
# associated with a factorization D = d_1*d_2 into negative fundamental
# discriminants.  Hence L_F(s,epsilon)=L(s,chi_{d_1})L(s,chi_{d_2}).
quadratic_genus_data = {
    12: (-3, -4, QQ(1) / 9, QQ(36)),
    21: (-3, -7, QQ(32) / 63, QQ(63) / 8),
    24: (-3, -8, QQ(2) / 3, QQ(6)),
    28: (-4, -7, QQ(8) / 7, QQ(7) / 2),
    33: (-3, -11, QQ(4) / 3, QQ(3)),
    44: (-4, -11, QQ(3), QQ(4) / 3),
}

print("\nQuadratic class groups and exact genus-character values")
for discriminant, (d1, d2, expected_L, expected_alpha) in (
    quadratic_genus_data.items()
):
    F = quadratic_by_discriminant[ZZ(discriminant)]
    assert F.class_group().order() == 1
    assert F.narrow_class_group().order() == 2

    chi1 = kronecker_character(d1)
    chi2 = kronecker_character(d2)
    L1 = exact_L_one_minus_k(chi1, 3)
    L2 = exact_L_one_minus_k(chi2, 3)
    L_value = L1 * L2
    alpha = QQ(4) / L_value

    print(
        "D =",
        discriminant,
        "(d1,d2) =",
        (d1, d2),
        "factors =",
        (L1, L2),
        "L_F(-2,epsilon) =",
        L_value,
        "alpha =",
        alpha,
    )

    assert d1 * d2 == discriminant
    assert L_value == expected_L
    assert alpha == expected_alpha


# Integrality alone eliminates D=21,28,44.  The ramified-dyadic lower bound
# |alpha| >= 2(2^3-1)=14 eliminates D=24.  D=12 is eliminated only after the
# exact character signs are retained; D=33 is eliminated by the split-prime
# coefficient argument in the accompanying proof.
assert quadratic_genus_data[21][3] not in ZZ
assert quadratic_genus_data[28][3] not in ZZ
assert quadratic_genus_data[44][3] not in ZZ
assert abs(quadratic_genus_data[24][3]) < 2 * (2**3 - 1)

# In D=12, the unique prime p above 2 is ordinary-principal, generated by
# 1+sqrt(3), but is not narrow-principal because that generator has norm -2.
F12 = quadratic_by_discriminant[12]
b12 = F12.gen()
if F12.defining_polynomial() == x**2 - 3:
    negative_norm_generator = 1 + b12
else:
    # This branch is kept for robustness if the defining polynomial above is
    # changed to an isomorphic model.
    negative_norm_generator = None
assert negative_norm_generator is not None
assert negative_norm_generator.norm() == -2
prime_above_two_12 = F12.ideal(negative_norm_generator)
assert prime_above_two_12.norm() == 2
assert F12.ideal(2) == prime_above_two_12**2

alpha12 = quadratic_genus_data[12][3]
for ell in range(1, 4):
    eta_at_p = (-1) ** ell
    omega_at_p = -1
    local_rhs = eta_at_p * 2 ** (ell - 1) * (1 - omega_at_p * 2**3)
    assert local_rhs != alpha12

# For ell >= 4 the magnitude is >= 72 > 36.
assert 9 * 2**3 > alpha12
# The magnitude equation can only suggest ell=3, and parity then gives the
# wrong sign: -36 rather than +36.
assert 9 * 2 ** (3 - 1) == alpha12
assert (-1) ** 3 * 9 * 2 ** (3 - 1) == -alpha12


# In D=33, verify that 2 splits and each prime has norm 2.  The proof uses
# multiplicativity and the fact that the product of their inverse narrow
# classes is [(2)]^{-1}=1.
F33 = quadratic_by_discriminant[33]
factor_two_33 = F33.ideal(2).factor()
assert len(factor_two_33) == 2
assert all(exponent == 1 for _, exponent in factor_two_33)
assert all(prime.norm() == 2 for prime, _ in factor_two_33)

print("\nAll exact E3 Eisenstein--cuspidal finite checks passed.")
