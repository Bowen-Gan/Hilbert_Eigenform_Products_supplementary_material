#!/usr/bin/env python3
"""Exact low-weight arithmetic for the D=12 weight-one identity.

Standard library only. The ring decompositions, adelic component matching,
one-cusp statements, and nonzero generator constants are published or
theoretical inputs. We count coefficients of those Hilbert series and
evaluate the two primitive Dirichlet L(0)-values with exact Fractions.
No modular form, global analytic statement, or ring theorem is computed.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def monomials(weights, degree):
    """Exponent tuples of weighted degree 'degree' in a polynomial ring."""
    if not weights:
        return [()] if degree == 0 else []
    first, *tail = weights
    require(first > 0, "Generator weights must be positive.")
    return [(exponent,) + suffix
            for exponent in range(degree // first + 1)
            for suffix in monomials(tail, degree - first * exponent)]


def hilbert_coefficients(weights, numerator_degrees, bound):
    """Truncated formal expansion, independent of monomial enumeration."""
    coefficients = [int(k == 0) for k in range(bound + 1)]
    for weight in weights:
        require(weight > 0, "Denominator weights must be positive.")
        for degree in range(weight, bound + 1):
            coefficients[degree] += coefficients[degree - weight]
    return [sum(coefficients[degree - shift]
                for shift in numerator_degrees if shift <= degree)
            for degree in range(bound + 1)]


def dirichlet_l_zero(table):
    """L(0,chi)=-sum_a chi(a) B_1(a/q), B_1(x)=x-1/2."""
    modulus = len(table)
    return -sum(Fraction(table[a % modulus]) *
                (Fraction(a, modulus) - Fraction(1, 2))
                for a in range(1, modulus + 1))


def quadratic_norm(a, b):
    return a * a - 3 * b * b


def verify():
    bound = 12
    up = hilbert_coefficients((2, 3, 4), (0, 11), bound)
    mixed_sl2 = hilbert_coefficients((1, 3, 4), (0, 10), bound)
    mixed_symmetric = hilbert_coefficients((1, 3, 4), (0,), bound)
    mixed_unit_invariant = hilbert_coefficients((1, 4, 6), (0,), bound)
    require((up[4], up[5]) == (2, 1), "Upper-upper dimensions failed.")
    require((mixed_sl2[4], mixed_sl2[5]) == (3, 3),
            "Mixed SL2 dimensions failed.")
    require(all(mixed_sl2[k] == mixed_symmetric[k] for k in (4, 5)),
            "Galois symmetry cannot be inferred in the target weights.")
    require((mixed_unit_invariant[4], mixed_unit_invariant[5]) == (2, 2),
            "Mixed unit-invariant dimensions failed.")

    bases = {}
    for component, weights, series in (
            ("upper_upper", (2, 3, 4), up),
            ("mixed_symmetric", (1, 3, 4), mixed_symmetric),
            ("mixed_unit_invariant", (1, 4, 6), mixed_unit_invariant)):
        bases[component] = {}
        for degree in (4, 5):
            basis = monomials(weights, degree)
            require(len(basis) == series[degree],
                    "Monomial and series counts disagree.")
            bases[component][str(degree)] = [list(item) for item in basis]

    # Exact inequalities/signs behind the manuscript's component matching.
    require(12 < 16, "The quadratic Minkowski bound is not below two.")
    require(2 not in {a * a % 3 for a in range(3)},
            "The negative Pell equation is not excluded modulo three.")
    require(quadratic_norm(2, 1) == 1, "The displayed unit norm failed.")
    require(quadratic_norm(0, 2) == -12,
            "The principal component's different has wrong norm sign.")
    require(quadratic_norm(6, 0) == 36,
            "The other component's shifted ideal has wrong norm sign.")

    # The primitive quadratic characters of discriminants -4 and -3.
    chi_minus_four = (0, 1, 0, -1)
    chi_minus_three = (0, 1, -1)
    for table in (chi_minus_four, chi_minus_three):
        modulus = len(table)
        for a in range(modulus):
            for b in range(modulus):
                require(table[(a * b) % modulus] == table[a] * table[b],
                        "The Dirichlet-character table is not multiplicative.")
        require(table[-1] == -1, "The Dirichlet character must be odd.")
    l_minus_four = dirichlet_l_zero(chi_minus_four)
    l_minus_three = dirichlet_l_zero(chi_minus_three)
    require(l_minus_four == Fraction(1, 2), "L(0,chi_-4) failed.")
    require(l_minus_three == Fraction(1, 3), "L(0,chi_-3) failed.")
    l_epsilon = l_minus_four * l_minus_three
    constants = {"principal_mixed": (1 + 1) * l_epsilon / 4,
                 "nonprincipal_upper_upper": (1 - 1) * l_epsilon / 4}
    require(constants == {"principal_mixed": Fraction(1, 12),
                          "nonprincipal_upper_upper": Fraction(0)},
            "The weight-one constant vector failed.")
    product_scalar = 1 / constants["principal_mixed"]
    require(product_scalar == 12, "The product normalization must be twelve.")

    dimensions = {}
    for degree in (4, 5):
        principal_modular = mixed_unit_invariant[degree]
        other_modular = up[degree]
        dimensions[str(degree)] = {
            "central_character": "1" if degree % 2 == 0 else "epsilon",
            "principal_mixed_modular": principal_modular,
            "principal_mixed_constant_term_rank": 1,
            "principal_mixed_cuspidal": principal_modular - 1,
            "nonprincipal_upper_upper_modular": other_modular,
            "nonprincipal_upper_upper_constant_term_rank": 1,
            "nonprincipal_upper_upper_cuspidal": other_modular - 1,
            "full_adelic_modular": principal_modular + other_modular,
            "full_adelic_eisenstein": 2,
            "full_adelic_cuspidal": principal_modular + other_modular - 2,
        }
    require(dimensions["4"]["full_adelic_cuspidal"] == 2,
            "The full S4 dimension failed.")
    require(dimensions["5"]["full_adelic_cuspidal"] == 1,
            "The full S5 dimension failed.")

    return {
        "status": "PASS",
        "verification_scope": "Exact Hilbert-series coefficient counts and rational Dirichlet L(0)-values, conditional on the stated mathematical inputs.",
        "python_version": sys.version.split()[0],
        "python_optimization": sys.flags.optimize,
        "field": "Q(sqrt(3))",
        "discriminant": 12,
        "sources": {
            "Aoki2006Dimensions": "Theorem 2.2, first and third identities, and Section 4.1; the final superscript ++ should read --.",
            "Dursthoff2016": "Section 3.2 and Theorem 8.1.",
            "Dirichlet_special_values": "Generalized Bernoulli formula at m=1.",
        },
        "hilbert_series_coefficients_through_weight_12": {
            "upper_upper_unit_invariant": up,
            "mixed_full_SL2": mixed_sl2,
            "mixed_Galois_symmetric": mixed_symmetric,
            "mixed_Galois_symmetric_and_unit_invariant": mixed_unit_invariant,
        },
        "monomial_exponents_in_displayed_generator_order": bases,
        "dimensions": dimensions,
        "dirichlet_character_tables": {"chi_minus_four": list(chi_minus_four),
                                       "chi_minus_three": list(chi_minus_three)},
        "special_values": {"L(0,chi_minus_four)": str(l_minus_four),
                           "L(0,chi_minus_three)": str(l_minus_three),
                           "L(0,epsilon)": str(l_epsilon)},
        "weight_one_constant_vector": {key: str(value)
                                       for key, value in constants.items()},
        "product_normalization_scalar": str(product_scalar),
        "not_verified": [
            "The published ring structures or the correction to Aoki's final superscript.",
            "Adelic/classical component identification and one-cusp constant-term rank.",
            "Artin factorization for the narrow class character epsilon.",
            "The divisor-involution proof that E_1 vanishes on every Fourier coefficient of the nonprincipal component.",
            "The existence, holomorphy, and Hecke-eigenform interpretation of h4,h5 and their product.",
        ],
    }


def main():
    if not __debug__ or sys.flags.optimize:
        print("ERROR: optimized Python is not permitted for this verifier.",
              file=sys.stderr)
        return 2
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the JSON certificate.")
    args = parser.parse_args()
    try:
        result = verify()
    except ArithmeticError as error:
        print("FAIL: " + str(error), file=sys.stderr)
        return 1
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
        print("PASS: exact D=12 low-weight counts and normalization; JSON written to "
              + str(args.output))
    else:
        print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
