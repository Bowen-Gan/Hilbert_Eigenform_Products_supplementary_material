"""Finite checks for the repaired degree-seven and degree-eight E2 boundary.

The complete corrected Bordeaux septic table is filtered here; completeness
of that table and the octic list remains an explicit published external input.
This script verifies
the defining-polynomial discriminants, dyadic and septic local factorizations,
finite Euler products through 199, and the inequalities used for the one
analytic survivor in each degree.
"""

from sage.all import *
from number_field_table_inputs import table_rows, SEPTIC_BOUND

proof.all(True)
if not __debug__:
    raise RuntimeError("verification assertions require Python without -O/-OO")
x = polygen(QQ)


def zeta_partial_at_two(field, bound=199):
    """Exact product of finite local Dedekind-zeta factors for p <= bound."""
    value = QQ(1)
    for p in prime_range(bound + 1):
        for prime_ideal, _ in field.ideal(p).factor():
            value *= 1 / (1 - QQ(prime_ideal.norm()) ** (-2))
    return value


def factor_degrees_mod_p(polynomial, p):
    return sorted([factor.degree() for factor, multiplicity in
                   polynomial.change_ring(GF(p)).factor()
                   for _ in range(multiplicity)])


# Degree eight: the complete external list below the analytic cutoff has
# discriminants 282300416, 309593125, and 324000000.
f8a = x**8 - 4*x**7 + 14*x**5 - 8*x**4 - 12*x**3 + 7*x**2 + 2*x - 1
f8b = x**8 - 4*x**7 - x**6 + 17*x**5 - 5*x**4 - 23*x**3 + 6*x**2 + 9*x - 1
f8c = x**8 - 7*x**6 + 14*x**4 - 8*x**2 + 1

F8a = NumberField(f8a, "a8")
F8b = NumberField(f8b, "b8")
F8c = NumberField(f8c, "c8")
assert [F8a.discriminant(), F8b.discriminant(), F8c.discriminant()] == [
    282300416, 309593125, 324000000
]
assert all(field.degree()==8 and field.signature()==(8,0)
           for field in [F8a,F8b,F8c])
assert [f8a.discriminant(), f8b.discriminant(), f8c.discriminant()] == [
    282300416, 309593125, 324000000
]
assert factor_degrees_mod_p(f8a, 2) == [4, 4]
assert factor_degrees_mod_p(f8c, 2) == [4, 4]
assert factor_degrees_mod_p(f8b, 2) == [8]
for field in [F8a,F8c]:
    dyadic=list(field.ideal(2).factor())
    assert len(dyadic)==1 and dyadic[0][1]==2 and dyadic[0][0].norm()==16
assert [(P.norm(),e) for P,e in F8b.ideal(2).factor()]==[(256,1)]

P8 = zeta_partial_at_two(F8b)
RBF = RealBallField(200)
assert RBF(P8) < RBF("1.010297")
assert RBF(P8) * (RBF(200) / 199) ** 8 < RBF("1.051634")
assert (4 * RBF.pi()**2) ** 8 / RBF(309593125) ** (RBF(3) / 2) > RBF("1.08315")


# Degree seven: after three fields are removed by the exact local data in
# Driver--Jones' table, the sole analytic survivor has discriminant 25367689.
f7 = x**7 - x**6 - 6*x**5 + 4*x**4 + 9*x**3 - 4*x**2 - 3*x + 1
F7 = NumberField(f7, "a7")
assert F7.discriminant() == 25367689
assert f7.discriminant() == 25367689
# Retain every source row under the analytic cutoff rounded upwards. The
# bundled source has complete range D<150000000 and preserves provenance.
septic_rows = table_rows(7, SEPTIC_BOUND)
assert [row["discriminant"] for row in septic_rows] == [20134393,25164057,25367689,28118369]
local_tests = {20134393:(7,[1,6]),25164057:(3,[1,1,5]),28118369:(7,[1,6])}
for source_row in septic_rows:
    D=ZZ(source_row["discriminant"])
    poly=PolynomialRing(QQ,"x")(source_row["polynomial"])
    F=NumberField(poly,"s")
    assert F.signature()==(7,0)
    assert F.discriminant()==D
    print("SEPTIC_FIELD",D,poly)
    if ZZ(D) in local_tests:
        p0,expected=local_tests[ZZ(D)]
        fac=list(F.ideal(p0).factor())
        degrees=sorted([ZZ(P.residue_class_degree()) for P,e in fac for j in range(e)])
        assert degrees==expected
        assert any(P.norm()**2<=2**7 and P**2!=F.ideal(2) for P,e in fac)
        print("SEPTIC_LOCAL",D,p0,degrees)
P7 = zeta_partial_at_two(F7)
assert RBF(P7) < RBF("1.022705")
assert RBF(P7) * (RBF(200) / 199) ** 7 < RBF("1.059227")
assert (4 * RBF.pi()**2) ** 7 / RBF(25367689) ** (RBF(3) / 2) > RBF("1.16975")

print("P7_EXACT =", P7)
print("P8_EXACT =", P8)
print("P7 =", RBF(P7))
print("P8 =", RBF(P8))
print("UPPER7 =", RBF(P7)*(RBF(200)/199)**7)
print("REQUIRED7 =", (4*RBF.pi()**2)**7 / RBF(25367689)**(RBF(3)/2))
print("UPPER8 =", RBF(P8)*(RBF(200)/199)**8)
print("REQUIRED8 =", (4*RBF.pi()**2)**8 / RBF(309593125)**(RBF(3)/2))
print("E2 high-degree finite checks passed")
