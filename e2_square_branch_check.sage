"""Exact quartic check for the square-ramified dyadic branch at weight two.

This is a finite certificate for the inertness proposition in Case 11.
It enumerates every totally real quartic field of discriminant below 4446,
selects precisely those for which (2)=p^2 with N(p)=4, and verifies that
each survivor has ordinary class number one.  It does not formalize the
analytic or Hilbert-modular-form arguments surrounding the enumeration.
"""

from sage.all import *
from sage.rings.number_field.totallyreal_rel import enumerate_totallyreal_fields_all

proof.all(True)
x = polygen(QQ)

fields = enumerate_totallyreal_fields_all(4, 4446, return_pari_objects=False)
square_rows = []

for discriminant, polynomial in fields:
    F = NumberField(polynomial, "a")
    factorization = list(F.ideal(2).factor())
    if (
        len(factorization) == 1
        and factorization[0][1] == 2
        and factorization[0][0].norm() == 4
    ):
        square_rows.append(
            (ZZ(discriminant), polynomial, ZZ(F.class_group().order()))
        )

expected_discriminants = [1600, 2000, 2624, 3600, 4400]
assert [row[0] for row in square_rows] == expected_discriminants
assert all(row[2] == 1 for row in square_rows)

print("quartic fields enumerated:", len(fields))
for discriminant, polynomial, class_number in square_rows:
    print(discriminant, polynomial, "h =", class_number)

# The rational comparison used after the enumeration is strict.
RBF = RealBallField(200)
assert RBF(4000) / RBF.pi() ** 8 > RBF(4) / 15

print("E2 square-ramified quartic certificate: all checks passed")
