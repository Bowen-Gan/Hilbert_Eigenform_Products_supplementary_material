#!/usr/bin/env python3
"""Independently verify the D=1125 fixed-character quaternion certificate.

Only Python's standard library is used. All field/quaternion arithmetic is
rebuilt from rational coefficients, rather than trusting Sage's recorded
character sum or central-character conclusion. This verifies finite inputs;
the mass formula, norm-class surjectivity, quaternionic decomposition, and
Jacquet--Langlands transfer remain explicitly cited mathematical theorems.
"""

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def determinant(rows):
    rows = [list(row) for row in rows]
    value = Fraction(1)
    for j in range(len(rows)):
        pivot = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if pivot is None:
            return value * 0
        if pivot != j:
            rows[pivot], rows[j] = rows[j], rows[pivot]
            value = -value
        entry = rows[j][j]
        value *= entry
        for i in range(j + 1, len(rows)):
            ratio = rows[i][j] / entry
            rows[i] = [u - ratio*v for u, v in zip(rows[i], rows[j])]
    return value


class FieldElement:
    """Q[a]/(a^4-a^3-4a^2+4a+1), with exact rational arithmetic."""

    def __init__(self, value=0):
        if isinstance(value, FieldElement):
            self.coefficients = value.coefficients
        elif isinstance(value, (list, tuple)):
            require(len(value) <= 4, "Field element has too many coefficients")
            self.coefficients = tuple(map(Fraction, value)) + (Fraction(0),)*(4-len(value))
        else:
            self.coefficients = (Fraction(value), Fraction(0), Fraction(0), Fraction(0))

    def __hash__(self):
        return hash(self.coefficients)

    def __bool__(self):
        return any(self.coefficients)

    def __eq__(self, other):
        return self.coefficients == FieldElement(other).coefficients

    def __neg__(self):
        return FieldElement([-c for c in self.coefficients])

    def __add__(self, other):
        return FieldElement([a+b for a, b in zip(self.coefficients, FieldElement(other).coefficients)])

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-FieldElement(other))

    def __rsub__(self, other):
        return FieldElement(other) - self

    def __mul__(self, other):
        other = FieldElement(other)
        coefficients = [Fraction(0)]*7
        for i, u in enumerate(self.coefficients):
            for j, v in enumerate(other.coefficients):
                coefficients[i+j] += u*v
        for degree in range(6, 3, -1):
            # a^4 = -1 - 4a + 4a^2 + a^3.
            entry = coefficients[degree]
            for j, multiplier in enumerate((-1, -4, 4, 1)):
                coefficients[degree-4+j] += multiplier*entry
        return FieldElement(coefficients[:4])

    __rmul__ = __mul__

    def __truediv__(self, rational):
        rational = Fraction(rational)
        return FieldElement([c/rational for c in self.coefficients])

    def __pow__(self, exponent):
        require(exponent >= 0, "Only nonnegative field powers are needed")
        result = FieldElement(1)
        for _ in range(exponent):
            result *= self
        return result

    def norm(self):
        images = [(self * FieldElement([0]*j + [1])).coefficients for j in range(4)]
        return determinant([[images[j][i] for j in range(4)] for i in range(4)])

    def serialized(self):
        return [str(c) for c in self.coefficients]


class Quaternion:
    def __init__(self, coefficients):
        require(len(coefficients) == 4, "Quaternion requires four field coefficients")
        self.coefficients = tuple(map(FieldElement, coefficients))

    def __hash__(self):
        return hash(self.coefficients)

    def __eq__(self, other):
        return isinstance(other, Quaternion) and self.coefficients == other.coefficients

    def __neg__(self):
        return Quaternion([-c for c in self.coefficients])

    def __add__(self, other):
        return Quaternion([a+b for a, b in zip(self.coefficients, other.coefficients)])

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        a, b, c, d = self.coefficients
        e, f, g, h = other.coefficients
        return Quaternion([a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                           a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e])

    def trace(self):
        return 2*self.coefficients[0]

    def norm(self):
        return sum(c*c for c in self.coefficients)

    def serialized(self):
        return [c.serialized() for c in self.coefficients]


def verify(data_dir):
    certificate_path = Path(data_dir)/"quartic_1125_sufficient.json"
    data = json.loads(certificate_path.read_text())
    require(data["schema"] == "D1125-fixed-central-character-endpoint-v2",
            "Rerun the updated D1125 Sage certificate first")
    require(data["polynomial"] == [1, 4, -4, -1, 1], "Wrong field polynomial")
    require(data["proof_all"] is True, "Sage arithmetic must use proof mode")
    a = FieldElement([0, 1])
    tau = a**3 - 3*a + 1
    require(tau**2-tau-1 == 0, "Icosian parameter polynomial fails")
    require(FieldElement(data["tau"]) == tau, "Recorded tau is incorrect")
    identity = Quaternion([1, 0, 0, 0])
    i = Quaternion([0, 1, 0, 0])
    z = Quaternion([tau/2, (tau-1)/2, Fraction(1, 2), 0])
    basis = [identity, i, z, i*z]
    require([Quaternion(q) for q in data["maximal_order_basis"]] == basis,
            "Incorrect recorded maximal order basis")
    gram = [[(u*v).trace() for v in basis] for u in basis]
    # The matrix from the manuscript gives the determinant algebraically.
    expected_gram = [[2, 0, tau, 1-tau], [0, -2, 1-tau, -tau],
                     [tau, 1-tau, tau-1, -1], [1-tau, -tau, -1, -tau]]
    require(gram == [[FieldElement(c) for c in row] for row in expected_gram],
            "Trace Gram matrix is not the manuscript matrix")
    # Its field determinant is computed by the permutation formula, avoiding
    # field division and independently checking the unit discriminant.
    from itertools import permutations
    gram_determinant = FieldElement(0)
    for permutation in permutations(range(4)):
        inversions = sum(permutation[j] > permutation[k]
                         for j in range(4) for k in range(j+1, 4))
        term = FieldElement((-1)**inversions)
        for j in range(4):
            term *= gram[j][permutation[j]]
        gram_determinant += term
    require(gram_determinant == -1, "Trace discriminant is not a unit")
    require(gram == [[FieldElement(c) for c in row] for row in data["trace_gram"]],
            "Recorded trace Gram matrix is incorrect")
    # Verify every recorded basis product directly and its integral coordinates.
    products = data["multiplication_table"]
    require(len(products) == 4 and all(len(row) == 4 for row in products),
            "Incomplete order multiplication table")
    for j, u in enumerate(basis):
        for k, v in enumerate(basis):
            require(len(products[j][k]) == 4, "Incomplete product coordinates")
            reconstructed = Quaternion([0, 0, 0, 0])
            for coefficient, b in zip(products[j][k], basis):
                c = FieldElement(coefficient)
                # Recorded coefficients are in Z[tau], hence algebraic integers
                # even though a power-basis coordinate can be rational.
                require(all(x.denominator == 1 for x in c.coefficients),
                        "Order coordinates are not integral power-basis elements")
                reconstructed += Quaternion([c*x for x in b.coefficients])
            require(reconstructed == u*v, "Order multiplication table fails")

    group = {identity}
    frontier = [identity]
    while frontier:
        u = frontier.pop()
        for generator in (i, z):
            v = u*generator
            if v not in group:
                group.add(v)
                frontier.append(v)
            require(len(group) <= 120, "Generated group has more than 120 elements")
    require(len(group) == 120 and -identity in group, "Incorrect binary group order")
    require(all(g.norm() == 1 for g in group), "Generated element is not norm one")
    rows = data["weight_five_character_rows"]
    require(len(rows) == 120, "Character table must have 120 rows")
    recorded = set()
    distribution = Counter()
    trace_distribution = Counter()
    for row in rows:
        g = Quaternion(row["element"])
        require(g in group and g not in recorded, "Character table repeats or invents an element")
        recorded.add(g)
        trace = g.trace()
        require(trace == FieldElement(row["trace"]), "Incorrect reduced trace")
        character = (trace**3-2*trace).norm()
        require(character.denominator == 1, "Nonintegral character")
        require(character == row["weight_five_character"], "Incorrect character value")
        distribution[int(character)] += 1
        trace_distribution[tuple(trace.serialized())] += 1
    require(recorded == group and distribution == {256: 2, 1: 88, 0: 30},
            "Wrong weight-five character distribution")
    character_sum = sum(value*count for value, count in distribution.items())
    require(character_sum == data["weight_five_character_sum"] == 600,
            "Incorrect character sum")
    require(data["principal_ideal_class_weight_five_invariants"] == character_sum//120 == 5,
            "Incorrect invariant dimension")
    require(data["weight_five_character_distribution"] ==
            [{"value":value,"multiplicity":count} for value,count in sorted(distribution.items())],
            "Incorrect recorded character distribution")
    require(data["reduced_trace_distribution"] ==
            [{"trace":list(t),"multiplicity":m} for t,m in sorted(trace_distribution.items())],
            "Incorrect recorded reduced trace distribution")
    units = [FieldElement(u) for u in data["fundamental_units"]]
    require(len(units) == 3 and all(u.norm() == 1 for u in units), "Base-unit scalar action fails")
    signatures = data["unit_signature_image"]
    from itertools import product
    even_signs = sorted(list(s) for s in product((-1, 1), repeat=4)
                        if s[0]*s[1]*s[2]*s[3] == 1)
    require(sorted(signatures) == even_signs and data["signature_image_order"] == 8,
            "Recorded unit-signature image is not the even-sign subspace")
    central = data["fixed_central_character"]
    require(central["name"] == "epsilon" and central["infinity_type"] == [1]*4,
            "Wrong target central character")
    require(central["base_unit_action"] == 1, "Wrong base-unit action")
    expected_scalars = {"1": FieldElement(1), "-1": FieldElement(-1),
                        "2": FieldElement(2), "a+1": a+1}
    expected_scalars.update({"unit_%s"%(j+1):u for j,u in enumerate(units)})
    require(len(central["scalar_rows"]) == len(expected_scalars) and
            {row["name"] for row in central["scalar_rows"]} == set(expected_scalars),
            "Incomplete scalar action examples")
    central_rows = []
    for row in central["scalar_rows"]:
        b = FieldElement(row["element"])
        require(b == expected_scalars[row["name"]], "Incorrect scalar representative")
        norm = b.norm()
        require(norm != 0 and norm == Fraction(row["field_norm"]), "Wrong scalar norm")
        action = norm**3 / abs(norm)**3
        require(Fraction(row["algebraic_weight_five_scalar_action"]) == norm**3,
                "Wrong algebraic scalar action")
        require(Fraction(row["principal_ideal_norm"]) == abs(norm), "Wrong ideal norm")
        require(Fraction(row["normalized_inverse_translation_scalar"]) == action,
                "Wrong normalized central action")
        require(row["epsilon_principal_ideal_value"] == (1 if norm > 0 else -1) == action,
                "Wrong epsilon ideal value")
        central_rows.append({"name": row["name"], "field_norm": str(norm),
                             "normalized_scalar": str(action)})
    require({row["normalized_scalar"] for row in central_rows} == {"-1", "1"},
            "Both central sign cosets must be checked")
    require((a+1).norm() == -5, "Negative norm representative fails")
    require(data["ordinary_class_number"] == 1 and data["narrow_class_number"] == 2,
            "Incorrect class-number inputs")
    require(data["quaternion_ideal_class_number"] == 2, "Wrong quaternion class count")
    require(data["norm_one_subgroup_order"] == 120 and
            data["full_reduced_unit_group_order"] == 60 and
            data["universal_reduced_unit_group_bound"] == 60,
            "Wrong reduced unit-group bound or order")
    from math import gcd
    possible_orders = [m for m in range(1,129)
                       if sum(gcd(j,m) == 1 for j in range(1,m+1)) <= 8]
    require(data["possible_projective_element_orders"] == possible_orders and
            max(possible_orders) == 30, "Incomplete projective element-order candidates")
    require(data["local_euler_factors"] == [{"p":2,"factors":[[1,4]]},
                                           {"p":3,"factors":[[2,2]]}],
            "Wrong local splitting inputs")
    local_ratio = Fraction(27,85)*Fraction(256,405)
    require(local_ratio/Fraction(2**7*6**4) == Fraction(data["mass_upper_coefficient"]),
            "Mass bound coefficient does not match the local Euler factors")
    mass_square = Fraction(data["mass_upper_squared"])
    require(Fraction(data["mass_upper_coefficient"]) == Fraction(1, 826200),
            "Incorrect rational mass upper coefficient")
    require(mass_square == Fraction(1125**3, 826200**2), "Incorrect squared mass bound")
    require(mass_square < Fraction(1, 400) and 60**2*mass_square < 9,
            "Mass bound no longer forces fewer than three ideal classes")
    require(data["conclusions"] == {"dim_S2_trivial_character":0,
                                     "dim_S5_totally_odd_character_at_least":5},
            "Incorrect stated dimension conclusions")
    return {
        "schema": "D1125-independent-quaternion-central-character-v1",
        "result": "PASS", "certificate": str(certificate_path),
        "field_arithmetic": "exact rational polynomial quotient",
        "norm_one_group_order": 120, "reduced_group_order": 60,
        "trace_gram_determinant": "-1", "weight_five_character_sum": 600,
        "character_distribution": [{"value":v,"multiplicity":m} for v,m in sorted(distribution.items())],
        "trace_distribution": [{"trace":list(t),"multiplicity":m} for t,m in sorted(trace_distribution.items())],
        "central_action_rule": "Norm(b)^3/abs(Norm(b))^3 = epsilon((b))",
        "central_scalar_checks": central_rows, "fixed_central_character": "epsilon",
        "dim_S5_epsilon_lower_bound": 5,
        "theorem_dependencies": data["theorem_dependencies"],
        "scope": "Independent finite verification; class-number and unit-signature completeness are proof-mode Sage inputs; representation, mass and Jacquet--Langlands statements are cited theorems.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent/"rerun")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.data_dir)
    output = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(output)
        print("PASS: quaternion arithmetic and epsilon central action; dim_S5_epsilon_lower_bound=5")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
