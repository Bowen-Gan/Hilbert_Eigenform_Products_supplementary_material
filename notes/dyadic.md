# Dyadic coefficient certificates for source weights at least two

`e2_dyadic_bound.sage` enumerates every ordered totally positive decomposition
of 4 and checks the normalized coefficient inequality used in the manuscript.
`verify_dyadic_bounds.py` independently reconstructs the same arithmetic using
Python integers and rational intervals. The current outputs belong in `rerun/`
and `verification_case11/`; run status and logs describe the actual execution.

| Discriminant | Bound A for the absolute reciprocal constant | Ordered decompositions of 4 | Excluded integer source weights |
|---|---:|---:|---|
| 49 | 168 (exact alpha = -168) | 9 | At least 7 |
| 169 | 24 (exact alpha = -24) | 3 | At least 4 |
| 361 | 8 (exact alpha = -8) | 3 | At least 3 |
| 725 | 121 | 7 | At least 4 |
| 1125 | 60 (exact alpha = 60) | 15 | At least 4 |
| 5125 | 61/10 | 7 | At least 2 |
| 6125 | Exact alpha = 60/13 | Not required | Excluded by nonintegrality |

At discriminants 49, 169, 361, and 725 the narrow class number is one,
so compatible source weights are even. These comparisons leave 2, 4, 6
at discriminant 49, only 2 at discriminants 169, 361, and 725, and 2, 3
at discriminant 1125. The remaining cusp-space arguments are separate
mathematical inputs: `cubic_weight_two_genus.sage`,
`quartic_725_weight_two.sage`, and `quartic_1125_sufficient.sage` check
arithmetic for the small-space arguments. The discriminant-49 modular-form
inputs come from Borisov--Gunnells, as stated in the manuscript.
The manuscript treats source weight one by its endpoint pairing and twisting
argument; none of the finite searches here is presented as an unresolved case.

## Exact special values

For each abelian field, an exact Gaussian-period isomorphism identifies the
cyclic quotient character group. The product of primitive Dirichlet special
values gives

    zeta_49(-1) = -1/21,
    zeta_169(-1) = -1/3,
    zeta_361(-1) = -1,
    zeta_1125(-1) = 4/15,
    zeta_6125(-1) = 52/15.

Each factor is a finite generalized Bernoulli sum; no numerical reconstruction
of a rational special value is used. `verify_special_values.py` independently
checks the first four values using rational coordinate arithmetic in
Q(zeta_3) and Q(i). It also checks the weight-three and weight-four
Dirichlet tables used elsewhere in the manuscript.

## Rational Euler-product bounds

For discriminants 725 and 5125, let P be the product over all prime ideals
above rational primes at most 199. Exact ideal factorization computes P over
QQ, and P is strictly smaller than zeta_F(2). The functional equation, pi < 22/7,
and rational lower bounds for sqrt(D) give

    |alpha| < (4*(22/7)^2)^4 / (D*sqrt_lower(D)*P).

The square-root lower bounds are 2692582403/10^8 and 7158910531/10^8,
respectively; integer squaring certifies their directions. The products satisfy
P > 1.036232 and P > 1.102613, respectively. The resulting reciprocal bounds
are |alpha| < 120.468280 and |alpha| < 6.023835, so the manuscript's choices
A=121 and A=61/10 are valid. Exact algebraic comparisons at the selected source
weights give right sides below 138 and 149, respectively, while the left side
is 255. The recorded certificate uses these same A values.

## Exhaustive enumeration and integer weight monotonicity

Every listed polynomial has polynomial discriminant equal to the field
discriminant, so its power basis is integral and maximal. Every x with
0 << x << 4 satisfies Tr(x^2) < 16*[F:Q]. Exact short-vector enumeration
uses a finite integer box derived from the inverse trace Gram matrix:
c_j^2 < 16*[F:Q]*(G^-1)_{jj}. Filtering the box by the strict trace bound
includes every short vector with both signs; exact total-positivity tests
retain precisely all such x. The independent Python implementation
reconstructs its own Gram matrix and exhaustively checks the same proof bound.

After division by q^(ell-1), where q=2^[F:Q], each convolution term has
base sqrt(N(4-x))/q strictly between zero and one; the other two bases are
q^(-1/2) and q^(-1). Increasing the integer source weight therefore decreases
every positive term. A successful cutoff comparison excludes every subsequent
integer weight. Decimal intervals are display data only.

Run from the repository root:

    sage -python e2_dyadic_bound.sage --output rerun/e2_dyadic_output.json
    python3 verify_dyadic_bounds.py --data-dir rerun --output verification_case11/dyadic_verification.json
