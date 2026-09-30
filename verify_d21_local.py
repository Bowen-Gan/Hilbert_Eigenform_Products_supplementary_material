#!/usr/bin/env python3
"""Exact local arithmetic checks for the Case 11 quadratic endpoints.

This is a standard-library-only certificate for the arithmetic in the
discriminant-21 exclusion and the discriminant-40 indecomposability test.
It uses integer inequalities, never decimal embeddings.  Run, for example,

    python verify_d21_local.py --output verification.json

The infinite weight exclusions are checked by their integer base cases and
an explicit induction in steps of two; they are not inferred from a finite
weight search.  Assertions are not used, and optimized Python is rejected.

External mathematical dependencies, which this program does NOT verify:
  * the normalization alpha = 12 at D=21 and h(F)=1, h+(F)=2;
  * the full-level Hecke recurrence and multiplicativity identities;
  * the Fourier product formula, including its component normalization;
  * character parity and the identification of the nontrivial narrow class;
  * the Ramanujan bound |c((2),h)| <= 2**ell for ell >= 2;
  * the positive-exponent norm inequality used for the primes above 5.
The output records these dependencies separately.  In particular, successful
execution does not certify the automorphic or Rankin--Selberg theorems.

Integer-arithmetic source snippets are included in the output to make the
finite enumeration and the ideal relations independently inspectable.
"""

import argparse
import json
from math import gcd, isqrt
from pathlib import Path
import sys


def check(condition, description):
    if not condition:
        raise ValueError("Verification failed: " + description)


def mul21(x, y):
    # w**2 = w + 5 in O_F = Z[w].
    a, b = x
    c, d = y
    return (a*c + 5*b*d, a*d + b*c + b*d)


def norm21(x):
    a, b = x
    return a*a + a*b - 5*b*b


def positive21(x):
    a, b = x
    twice_rational_part = 2*a + b
    return twice_rational_part > 0 and twice_rational_part**2 > 21*b*b


def decompositions21(m):
    # If 0 < x_sigma < m at both embeddings, |b| sqrt(21) < m.
    # B is the largest integer allowed by 21*b*b < m*m.
    bound = isqrt((m*m - 1)//21)
    result = []
    for b in range(-bound, bound + 1):
        # 0 < 2*a+b < 2*m gives the complete interval for a.
        first_a = -((b - 1)//2)
        last_a = (2*m - 1 - b)//2
        for a in range(first_a, last_a + 1):
            x, y = (a, b), (m - a, -b)
            if positive21(x) and positive21(y):
                result.append((x, y))
    return sorted(result)


def multiplication_matrix21(x):
    a, b = x
    return [[a, 5*b], [b, a + b]]


def determinant2(matrix):
    return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]


def verify_d21():
    expected = {
        3: [((1, 0), (2, 0)), ((2, 0), (1, 0))],
        4: [((1, 0), (3, 0)), ((2, 0), (2, 0)), ((3, 0), (1, 0))],
        5: [((1, 0), (4, 0)), ((2, 0), (3, 0)), ((3, 0), (2, 0)),
            ((4, 0), (1, 0)), ((2, 1), (3, -1)), ((3, -1), (2, 1))],
    }
    decompositions = {}
    for m, pairs in expected.items():
        actual = decompositions21(m)
        check(actual == sorted(pairs), "complete D21 decomposition list at " + str(m))
        decompositions[str(m)] = {
            "ordered_count": len(actual),
            "pairs_in_basis_1_w": actual,
            "complete_b_bound": isqrt((m*m - 1)//21),
        }

    unit, inverse = (2, 1), (3, -1)
    check(norm21(unit) == norm21(inverse) == 1, "the two norm-one units")
    check(positive21(unit) and positive21(inverse), "total positivity of the units")
    check(mul21(unit, inverse) == (1, 0), "the units are mutual inverses")
    p3_generator = (1, 1)
    check(norm21(p3_generator) == -3, "prime-above-3 generator norm")
    check(mul21(p3_generator, p3_generator) == (6, 3), "p3 generator squared")
    check((6, 3) == tuple(3*x for x in unit), "p3 squared is 3 times a unit")
    check(all((r*r - r - 5) % 3 == (r + 1)**2 % 3 for r in range(3)),
          "repeated factor of the defining polynomial modulo 3")
    check(all((r*r - r - 5) % 2 != 0 for r in range(2)), "2 is inert")

    p5_generator, q5_generator = (0, 1), (1, -1)
    check(norm21(p5_generator) == norm21(q5_generator) == -5,
          "the two prime-above-5 generator norms")
    check(mul21(p5_generator, q5_generator) == (-5, 0), "p5*q5=(5)")
    check(tuple(x + y for x, y in zip(p5_generator, q5_generator)) == (1, 0),
          "coprimality of the two primes above 5")
    for x in (p3_generator, p5_generator, q5_generator):
        check(determinant2(multiplication_matrix21(x)) == norm21(x),
              "principal-ideal index equals absolute norm")

    # With a2=-5-2*(-1)**ell*3**(ell-2), Ramanujan allows |a2|<=2**ell.
    # Even ell=2 is checked separately; ell=4 starts step-two induction.
    check(5 + 2*3**0 > 2**2, "even-weight base ell=2")
    check(5 + 2*3**2 > 2**4, "even-weight base ell=4")
    # For t=3**(ell-2)>=9: A(ell+2)-4*A(ell)=10*t-15>0.
    check(10*9 - 15 > 0, "even-weight induction minimum")
    check(2*3**3 - 5 > 2**5, "odd-weight base ell=5")
    # For odd ell>=5: A(ell+2)-4*A(ell)=10*t+15>0.
    check(10 > 0 and 15 > 0, "odd-weight induction coefficients")

    ell = 3
    a2 = -5 - 2*(-1)**ell*3**(ell - 2)
    a3 = -1 - 5*4**(ell - 2) - 3*a2
    a4 = a2*a2 - 4**(ell - 1)
    check((a2, a3, a4) == (1, -24, -15), "remaining-weight cusp coefficients")
    e2, e3, e4 = sum([1, 4]), sum([1, 3, 9]), sum([1, 4, 16])
    check((e2, e3, e4) == (5, 13, 21), "Eisenstein divisor sums")
    target_a2 = a2 + 12
    target_a3 = a3 - 8*(-1)**ell*3**(ell - 1)
    target_a4 = target_a2*target_a2 - 4**(ell + 1)
    check(target_a3 - a3 == 12*(a2 + e2),
          "agreement of the supplied Hecke and Fourier formulas at (3)")
    check(target_a4 - a4 == 12*(a3 + e2*a2 + e3),
          "agreement of the supplied Hecke and Fourier formulas at (4)")
    cusp_coefficients = {1: 1, 2: a2, 3: a3, 4: a4}
    eisenstein_coefficients = {1: 1, 2: e2, 3: e3, 4: e4}
    convolution_rows = []
    for x, y in decompositions["5"]["pairs_in_basis_1_w"]:
        e_x = eisenstein_coefficients[x[0]] if x[1] == 0 else 1
        a_y = cusp_coefficients[y[0]] if y[1] == 0 else 1
        if x[1] != 0:
            check(norm21(x) == 1, "nonrational Eisenstein exponent is a unit")
        if y[1] != 0:
            check(norm21(y) == 1, "nonrational cusp exponent is a unit")
        convolution_rows.append({"Eisenstein_exponent": x, "cusp_exponent": y,
                                 "Eisenstein_coefficient": e_x, "cusp_coefficient": a_y,
                                 "contribution": e_x*a_y})
    summands = [row["contribution"] for row in convolution_rows]
    check(len(summands) == 6 and sum(summands) == a4 + e2*a3 + e3*a2 + e4 + 2,
          "one convolution contribution for every ordered decomposition")
    product_difference = 12*sum(summands)
    check(product_difference == -1188, "Fourier-product coefficient difference")
    # Under the stated class/component hypotheses, norms of summand ideals
    # are >=3. The norm inequality would require 5 >= (2*sqrt(3))**2=12.
    check(5 < 12, "prime-above-5 indecomposability integer comparison")
    multiplicative_difference = 0
    check(product_difference != multiplicative_difference, "nonzero local contradiction")

    return {
        "field": {"discriminant": 21, "basis": ["1", "w"],
                  "defining_relation": "w^2-w-5=0"},
        "decompositions": decompositions,
        "units": {"unit": unit, "inverse": inverse, "norms": [1, 1],
                  "product": mul21(unit, inverse), "totally_positive": True},
        "prime_witnesses": {
            "p3_generator": p3_generator, "p3_generator_norm": -3,
            "p3_generator_square": [6, 3], "p3_square_relation": "(1+w)^2=3(2+w)",
            "p5_generator": p5_generator, "q5_generator": q5_generator,
            "p5_q5_product": [-5, 0], "p5_q5_sum": [1, 0],
            "2_inert_by_no_root_mod_2": True,
        },
        "infinite_weight_exclusion": {
            "coefficient_formula": "a2=-5-2*(-1)^ell*3^(ell-2)",
            "Ramanujan_bound_input": "abs(a2)<=2^ell, ell>=2",
            "even_bases": [{"ell": 2, "abs_a2": 7, "bound": 4},
                           {"ell": 4, "abs_a2": 23, "bound": 16}],
            "even_step": "A(ell+2)-4*A(ell)=10*t-15>0 for t=3^(ell-2)>=9",
            "even_step_minimum": 10*9 - 15,
            "odd_base": {"ell": 5, "abs_a2": 49, "bound": 32},
            "odd_step": "A(ell+2)-4*A(ell)=10*t+15>0 for t=3^(ell-2)>0",
            "only_unexcluded_source_weight_at_least_two": 3,
        },
        "remaining_weight": {
            "ell": ell, "cusp_coefficients_at_2_3_4": [a2, a3, a4],
            "target_coefficients_at_2_3_4_from_supplied_Hecke_formulas":
                [target_a2, target_a3, target_a4],
            "differences_at_3_4_from_both_supplied_formulas":
                [target_a3 - a3, target_a4 - a4],
            "Eisenstein_coefficients_at_2_3_4": [e2, e3, e4],
            "ordered_convolution_contributions_at_5": convolution_rows,
            "alpha_input": 12, "Fourier_product_difference_at_5": product_difference,
            "multiplicativity_difference_at_5_conditional_on_inputs": multiplicative_difference,
            "prime_indecomposability_comparison": "5<12=(sqrt(3)+sqrt(3))^2",
        },
        "source_snippets": [
            "x=a+bw=(2a+b+b*sqrt(21))/2; total positivity is 2a+b>0 and (2a+b)^2>21b^2.",
            "For 0<x_sigma<m at both embeddings, 21b^2<m^2 and 0<2a+b<2m; these bound every integer pair enumerated.",
            "Norm(a+bw)=a^2+ab-5b^2 is the determinant of multiplication by a+bw on Z[w].",
            "(1+w)^2=3(2+w), (2+w)(3-w)=1, w(1-w)=-5, w+(1-w)=1.",
        ],
    }


def positive40(x):
    a, b = x
    return a > 0 and a*a > 10*b*b


def verify_d40():
    nu = (7, 2)
    check(positive40(nu) and 7*7 - 10*2*2 == 9, "D40 nu positive of norm 9")
    # Both x and nu-x positive imply integer 1<=a<=6. Then 10*b^2<a^2<=36,
    # giving |b|<=1; the following 18 cases therefore cover every summand.
    b_bound = isqrt((6*6 - 1)//10)
    check(b_bound == 1, "complete D40 b bound")
    cases, positive_pairs = [], []
    for a in range(1, 7):
        for b in range(-b_bound, b_bound + 1):
            x, y = (a, b), (7 - a, 2 - b)
            x_positive, y_positive = positive40(x), positive40(y)
            cases.append({"x": x, "nu_minus_x": y,
                          "x_totally_positive": x_positive,
                          "nu_minus_x_totally_positive": y_positive})
            if x_positive and y_positive:
                positive_pairs.append((x, y))
    check(len(cases) == 18 and not positive_pairs, "D40 complete indecomposability check")

    # p=(3,sqrt(10)-1); p^2 is spanned by 9, 3(sqrt(10)-1), (sqrt(10)-1)^2.
    square_generators = [(9, 0), (-3, 3), (11, -2)]
    check(determinant2([[3, -1], [0, 1]]) == 3, "D40 prime-ideal lattice index")
    minors = [square_generators[i][0]*square_generators[j][1]
              - square_generators[i][1]*square_generators[j][0]
              for i in range(3) for j in range(i + 1, 3)]
    lattice_index = gcd(*minors)
    check(lattice_index == 9, "D40 square-ideal lattice index")
    combination = [-1, 2, 2]
    generated_nu = tuple(sum(c*g[k] for c, g in zip(combination, square_generators))
                         for k in range(2))
    check(generated_nu == nu, "D40 nu lies in p squared")
    quotients = []
    for u, v in square_generators:
        rational_numerator, radical_numerator = 7*u - 20*v, 7*v - 2*u
        check(rational_numerator % 9 == radical_numerator % 9 == 0,
              "D40 each square-ideal generator is divisible by nu in Z[sqrt(10)]")
        quotients.append([rational_numerator//9, radical_numerator//9])
    check([r*r % 3 for r in range(3)] == [0, 1, 1], "D40 roots modulo 3")

    return {
        "field": {"discriminant": 40, "basis": ["1", "sqrt(10)"]},
        "nu": nu, "norm_nu": 9, "nu_totally_positive": True,
        "complete_a_range": [1, 6], "complete_b_range": [-1, 1],
        "enumerated_cases": cases, "positive_decompositions": positive_pairs,
        "nu_indecomposable": True,
        "prime_square_witness": {
            "prime_generators": [[3, 0], [-1, 1]],
            "prime_norm": 3, "square_generators": square_generators,
            "maximal_minors": minors, "square_lattice_index": lattice_index,
            "nu_linear_combination_coefficients": combination,
            "square_generator_quotients_by_nu": quotients,
            "ideal_relation": "(3,sqrt(10)-1)^2=(7+2sqrt(10))",
        },
        "source_snippets": [
            "For x=a+b*sqrt(10), total positivity is a>0 and a^2>10b^2.",
            "If x and nu-x are totally positive, 1<=a<=6 and 10b^2<a^2<=36, so b is -1,0,1.",
            "-9+2*3*(sqrt(10)-1)+2*(sqrt(10)-1)^2=7+2sqrt(10).",
            "The gcd of maximal 2x2 minors of the three square-ideal generators is 9, equal to Norm(nu).",
        ],
    }


def main():
    if sys.flags.optimize:
        raise SystemExit("Optimized Python (-O or -OO) is not accepted for this verifier.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", help="write the exact verification JSON to this path")
    args = parser.parse_args()
    result = {
        "schema": "case11-D21-D40-integer-local-verification-v1",
        "arithmetic": "Python standard-library integers; strict inequalities are exact",
        "D21": verify_d21(), "D40": verify_d40(),
        "external_dependencies_not_verified": [
            "D21 alpha=12, ordinary class number one and narrow class number two (separate rq_certificate.py).",
            "Full-level Hecke recurrences, multiplicativity, and Fourier product/component formulas.",
            "Character parity identifies the primes of negative generator norm with the nontrivial narrow class.",
            "Ramanujan bound for full-level parallel source weight ell>=2 (Blasius2006Ramanujan in manuscript).",
            "Positive-exponent norm inequality and same-component summand ideal-class identification.",
        ],
        "automorphic_theorem_verified": False,
        "all_integer_arithmetic_checks_passed": True,
    }
    serialized = json.dumps(result, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(serialized, encoding="utf-8")
        print("D21/D40 exact local arithmetic passed. Output: " + args.output)
    else:
        print(serialized, end="")


if __name__ == "__main__":
    main()
