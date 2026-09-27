# Case 11 dyadic calculation: completed and executed

The replacement `free_certificates/e2_dyadic_bound.sage` was executed with
SageMath 10.8.12 (passagemath) and all assertions passed. Its actual output is
`free_certificates/e2_dyadic_output.json`. It covers the exceptional discriminants
49, 1125, 5125, and 6125, and the additional lattice-search survivors 169,
361, and 725. The original uploads are unchanged.

| Discriminant | Certified input | Positive decompositions of 4 | Dyadic conclusion |
|---|---|---:|---|
| 49 | Exact alpha = -168 | 9 | Every source weight at least 7 is excluded |
| 169 | Exact alpha = -24 | 3 | Every source weight at least 4 is excluded |
| 361 | Exact alpha = -8 | 3 | Every source weight at least 3 is excluded |
| 725 | Absolute alpha < 120.08107414 | 7 | Every source weight at least 4 is excluded |
| 1125 | Exact alpha = 60 | 15 | Every source weight at least 4 is excluded |
| 5125 | Absolute alpha < 6.004472943 | 7 | Every source weight at least 2 is excluded |
| 6125 | Exact alpha = 60/13 | Not required | Alpha is not an algebraic integer |

For discriminant 49, the even-weight condition then leaves source weights
2, 4, and 6. The dyadic calculation by itself does not exclude weight 2.
The manuscript's phrase "the two remaining source weights" also uses
the absence of a weight-two cusp source. This separate modular-form input
is supplied by the Borisov--Gunnells citation, as explained in the main
manuscript additions: the three level-p weight-two cusp generators have
nontrivial g_7 eigenvalues zeta_7^3, zeta_7^6, zeta_7^5, hence the full-level
invariant space vanishes. It is not a separate output of this file.

For the additional discriminants 169, 361, and 725 the script also certifies
narrow class number one. Odd source weights are consequently unavailable.
Their dyadic inequalities leave only source weight 2. The separate free certificates cubic_weight_two_genus.sage and
quartic_725_weight_two.sage now prove S_2=0 for all three fields.
Thus these endpoints are closed without Magma.

## Exact special values

The script identifies the field as a Gaussian-period subfield of a
cyclotomic field, using an exact number-field isomorphism test. If chi
generates the quotient character group, it evaluates every primitive
character underlying chi^j by the finite formula

    L(-1, chi) = -1/2 sum(a=1..f) chi(a) (a^2/f - a + f/6).

The product over the quotient character group gives the Dedekind special
value. The computations give

    zeta_49(-1)   = -1/21,
    zeta_169(-1)  = -1/3,
    zeta_361(-1)  = -1,
    zeta_1125(-1) = 4/15,
    zeta_6125(-1) = 52/15.

In particular, the quartic factors of conductors 15 and 35 are respectively
the conjugate pairs -2 +/- 2i and 2 +/- 10i; their quadratic factor is -2/5
and the rational factor is -1/12. All field-identification data, kernel
residues, primitive conductors and character values are in the JSON output.
No numerical value is rounded to recover a rational special value.

Discriminants 169 and 361 are identified as the cubic Gaussian-period
fields of conductors 13 and 19 by the same exact method.

## Discriminant 5125: no exact special value needed

The functional equation gives

    |alpha| = (4*pi^2)^4 / (5125^(3/2) * zeta_F(2)).

The finite product over every prime ideal above rational primes at most
199 is a strict lower bound for zeta_F(2). This product is computed over QQ.
An outward-rounded 256-bit MPFI upper bound for pi and exact algebraic
arithmetic give the rational bound

    |alpha| < 6004472943/1000000000.

Using this upper bound in the normalized dyadic inequality at source
weight 2 gives

    RHS < 144.859588354 < 255 = 16^2 - 1.

Consequently the exact value alpha = 6 is unnecessary for this exclusion.
The old `e2_exact_zeta_values.m` is not a dependency of the replacement
certificate and can be omitted from the final supplement.

The same finite-Euler-product method for discriminant 725 gives
|alpha| < 6004053707/50000000 = 120.08107414. At source weight 4 its
normalized right side is below 136.148985297 < 255. At weight 3 the
bound is above 255, so the certified cutoff is 4.

## Exhaustive enumeration and the infinite weight range

For each relevant field the script computes an integral basis and its
positive trace Gram matrix. Every algebraic integer x with 0 << x << 4
satisfies Tr(x^2) < 16*[F:Q]. Exact short-vector enumeration below that
bound (including both signs), followed by exact total-positivity tests,
therefore finds all such x. The output records every coordinate vector,
the two algebraic integers, their norms, ideal-factor norms, divisor
counts and sigma_1 values.

After division by q^(ell-1), with q=2^[F:Q], each convolution summand has
base sqrt(N(4-x))/q < 1. The other terms also have bases below one.
Thus once the normalized right side is strictly smaller than q^2-1,
it remains so for all larger weights. All comparisons use Sage's exact
algebraic real field AA; decimal intervals are display data only.

## Minimal manuscript adjustment

The present conclusion at discriminant 5125 remains unchanged. Add a
short explanation such as:

> At discriminant 5125, a finite Euler product through 199 and the
> functional equation give |alpha| < 6.004472943. Substitution into the
> normalized dyadic inequality at source weight two gives a right side
> smaller than 144.86, whereas the left side is 255. Monotonicity excludes
> every larger source weight.

For the other two dyadic thresholds the actual comparisons are:

* D=49: at ell=6 the right side is 64.3753898187... > 63;
  at ell=7 it is exactly 33.9071044921875 < 63.
* D=1125: at ell=3 the right side is 286.40625 > 255;
  at ell=4 it is exactly 81.884765625 < 255.

Run the certificate with:

    sage -python e2_dyadic_bound.sage --output e2_dyadic_output.json

The equivalent `python` command works when `sage.all` is available in that
Python environment, as in the execution performed here.
