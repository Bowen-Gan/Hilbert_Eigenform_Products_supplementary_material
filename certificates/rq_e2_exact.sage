"""Independent Sage verification of the exact quadratic certificate.

Run from any directory with `sage -python rq_e2_exact.sage`, or with
`python rq_e2_exact.sage` in a Python environment providing sage.all.
Only arithmetic used in Case 11 is checked here.  D=21 is completed by the pure coefficient proof in notes/d21_local_free.tex.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rq_certificate as certificate
from sage.all import *
from sage.version import version as sage_version

proof.all(True)
print("SAGE_VERSION", sage_version)
reduction = certificate.reduction_certificate()
local = certificate.local_certificate()


def field_from_discriminant(D):
    x = polygen(QQ)
    polynomial = (x**2-x+(1-D)//4) if D % 4 == 1 else (x**2-D//4)
    return NumberField(polynomial, "w")


for row in reduction["rows"]:
    D = row["D"]
    F = field_from_discriminant(D)
    w = F.gen()
    assert F.discriminant() == D
    assert F.class_number(proof=True) == row["class_number"]
    assert F.narrow_class_group().order() == row["narrow_class_number"]
    value = kronecker_character(D).bernoulli(2)/24
    assert value == QQ(certificate.zeta_minus_one(D))
    for relation in row["relations"]:
        p, r = relation["P"]
        ideal = F.ideal(p, w-r)
        assert ideal.norm() == p
        kind = relation["kind"]
        if kind in ("base_squared", "prime_times_base"):
            q, s = row["base_ideal"]
            ideal *= F.ideal(q, w-s)
        a, b = relation["generator"]
        assert F.ideal(a+b*w) == ideal
        assert (a+b*w).norm() == relation["norm"]
    if row["class_number"] == 2:
        nontrivial = (kronecker_character(5).bernoulli(2)*
                      kronecker_character(D//5).bernoulli(2)/4)
        assert 4/nontrivial == QQ(row["alpha"][1])
    if row["negative_norm_unit"]:
        a, b = row["unit_certificate"]["unit"]
        assert (a+b*w).norm() == -1
    print("FIELD", D, "CLASS", row["class_number"], "NARROW",
          row["narrow_class_number"], "ZETA_MINUS_ONE", value,
          "RELATIONS", len(row["relations"]))

assert reduction["survivors"] == [12,21,24,28,69,77]
print("REDUCTION_SURVIVORS", reduction["survivors"])

for D in reduction["survivors"]:
    F = field_from_discriminant(D)
    factor2 = list(F.ideal(2).factor())
    P, e = factor2[0]
    assert len(factor2) == 1
    if D in [12,24,28]:
        assert e == 2 and P.norm() == 2
    else:
        assert e == 1 and P.norm() == 4
    print("DYADIC", D, factor2)

R = PolynomialRing(QQ, names=("h2", "h3"))
h2, h3 = R.gens()
for row in local["inert"]:
    D = row["D"]
    F = field_from_discriminant(D)
    for m in [1,2,3]:
        factors = sorted((int(P.norm()), int(e)) for P,e in F.ideal(m).factor())
        assert factors == sorted(row["ideal_data"][m]["factor_norms_exponents"])
    alpha = QQ(row["alpha"])
    f3 = row["f3"]
    for ell in range(2,7):
        constant = alpha-f3-(15/alpha)*4**(ell-1)
        from_hecke = (h2+alpha)**2-4**(ell+1)
        from_convolution = h2**2-4**(ell-1)+alpha*(h3+5*h2+f3)
        assert from_hecke-from_convolution == -alpha*(3*h2+h3-constant)
    print("INERT", D, "F3", f3, "RELATION", row["coefficient_relation"],
          "BOUND_AT_T", row["terminal_t"], row["normalized_upper_bound"], "< 15")
    if D in [69,77]:
        print("WEIGHT_TWO", D, "CONSTANT", row["weight_two_constant"],
              "RAMANUJAN_BOUND", row["ramanujan_bound"])

F = field_from_discriminant(21)
w = F.gen()
P3, e = list(F.ideal(3).factor())[0]
assert e == 2 and P3.norm() == 3
assert F.ideal(1+w) == P3 and (1+w).norm() == -3
print("D21_PRIME3", P3, "GENERATOR", 1+w, "NORM", (1+w).norm())
print("PASS: exact reduction, independent Sage arithmetic, and local identities")
