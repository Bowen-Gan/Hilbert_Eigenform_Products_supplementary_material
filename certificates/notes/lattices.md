# Case 11 finite enumeration and component-lattice calculations

These calculations have actually been run in the installed Sage-compatible
`passagemath-standard 10.8.12` environment.  The `.sage` files deliberately use
ordinary Python syntax and can be run with either Sage or that Python
environment.  Input uploads have not been changed.

## Deliverables and execution

```sh
python free_certificates/e2_candidate_enumeration.sage --output free_certificates/e2_candidates.json
python free_certificates/e2_indecomposable_search.sage --input free_certificates/e2_candidates.json --output free_certificates/e2_lattices.json
python free_certificates/e2_allclass_balanced.sage --output free_certificates/d3969_certificate.json
```

The search exits with code 2 when it has visited all fields but has not
obtained a lattice certificate for every field.  That is an intentional
unresolved-search indication, not a failed algebraic assertion.  The
`search_complete` flag records completion of the prescribed search;
`all_nonexception_candidates_certified` records its strictly narrower
mathematical outcome.  External dyadic or dimension certificates must be
combined separately with this output, without relabelling their cases as
lattice exclusions.

## Complete field input

The manuscript's general even-ordinary-class-number argument is applied
before enumeration.  Under the stated necessary condition `n*h <= 14`, odd
`h>1` can occur only in degrees 3 and 4.  Therefore the complete remaining
search uses these inclusive, rigorously rounded-up discriminant bounds:

| Degree | Enumeration endpoint | Reason |
|---|---:|---|
| 3 | 4218 | General reciprocal-constant bound |
| 4 | 68133 | General reciprocal-constant bound |
| 5 | 209508 | Odd class number must be 1 |
| 6 | 2429064 | Odd class number must be 1 |

The endpoints are checked with outward 256-bit interval arithmetic against
`D^3 < (2*pi^4/3)^(2*n)` and `D^3 < (4*pi^2)^(2*n)`, respectively.
Fields of class number 1 in degrees 3 and 4 are subsequently subjected to
the sharper bound as well.  Degrees at least 7 are handled by the separate
high-degree argument, not this enumeration.

Sage's `enumerate_totallyreal_fields_all` enumerates primitive and
imprimitive fields.  Its official documentation is
<https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/totallyreal_rel.html>.
Every defining polynomial and every filtering decision is saved.  Fields
are never indexed solely by discriminant.  In particular the two quartic
fields of discriminant 64525 are both present, as are the two cubic fields
of discriminant 3969.

The completed run enumerates 771 fields: 143 cubics, 551 quartics, 37
quintics, and 40 sextics.  The stage counts are:

| Recorded outcome | Number |
|---|---:|
| Ordinary class number even | 15 |
| Excluded by the sharper class-number-one discriminant bound | 520 |
| `(2)` not inert | 135 |
| Small-prime norm criterion | 43 |
| Reciprocal constant cannot be an integer | 37 |
| Requires a lattice or exceptional arithmetic calculation | 21 |

The only surviving odd ordinary class number greater than one is 3, for
the cubic field `x^3 - 21*x - 35` of discriminant 3969.  This is now an
enumerated result, not an assumed list.

## Finite Euler filter

The old hardcoded list omitted fields without a complete filtering record.
The replacement computes all factor degrees through 199 and the exact
rational finite Euler product `P`.  For degree `n`,

`P < zeta_F(2) < P*(200/199)^n`.

The upper bound follows from the telescope
`product_(m>199) (1-m^-2)^-1 = 200/199`, after bounding the
rational-prime factor by `(1-p^-2)^-n`.  If `h=1`, then

`|alpha| = (4*pi^2)^n / (D^(3/2)*zeta_F(2))`.

The resulting certified interval is intersected with the integers.  An
empty intersection proves the required nonintegrality.  The transcript
contains all local factor degrees, the exact rational product and tail,
and outward numerical endpoints.  No candidate rational special value is
guessed or reconstructed.  This excludes, among others, discriminants
1369 and 6125.

## Exact all-component algorithm

For each integral ideal `A=P` or `P^2`, let `m` be the least positive
integer in `A`, and use `t=(m)*A^-1`.  This is an integral component lattice
with `A=(m)*t^-1`.  Scaling both exponent and component lattice by a
totally positive factor preserves decompositions, so no choice of
fundamental units and no narrow-principality assumption is needed.

An exactly checked unimodular LLL basis change reduces the trace Gram
matrix.  Trace-dual coordinates give a finite box containing every
`0 << x << m`.  Bounds use outward real intervals and rational endpoints;
signs use intervals when decisive and exact algebraic embeddings otherwise.
Every point in the containing box is tested.  A failed bounded search
records actual decomposition witnesses and is not treated as an exclusion.

The complete search through rational prime 43 finds these certificates:

| Degree | Discriminant | Norm of certificate prime |
|---|---:|---:|
| 3 | 81 | 3 |
| 3 | 257 | 3 |
| 3 | 321 | 3 |
| 3 | 697 | 5 |
| 3 | 1257 | 3 |
| 3 | 1489 | 7 |
| 3 | 3969 | 3 |
| 4 | 2525 | 5 |
| 4 | 4205 | 5 |
| 4 | 4525 | 5 |
| 4 | 8069 | 5 |
| 4 | 16317 | 5 |
| 5 | 14641 | 11 |
| 5 | 38569 | 7 |
| 6 | 371293 | 13 |

For each row the machine-readable output contains the defining polynomial,
prime ideal, component bases, trace Gram matrices, coordinate bounds, and
exhaustion counts for both `P` and `P^2`.

The preassigned external arithmetic cases are 49, 1125 and 5125; 6125
has already been excluded by the Euler filter.  **The additional fields
169, 361 and 725 have no such certificate for rational primes through 43.**
Their decomposition witnesses are retained.  Their dyadic/space treatment
must be included in the combined proof; the uploaded endpoint-space script
already mentions these fields.

## Discriminant 3969

The standalone calculation verifies the field discriminant, both class
groups `C3`, `P=(3,a+1)`, `Norm(P)=3`, and `(3)=P^3`.
For `P`, the component is `P^2`; for `P^2`, it is `P`.  The exponent is
3 in each.  In the recorded reduced bases both exact containing boxes
are `[0,1] x [-1,1] x [-1,1]`, with 18 points and no decomposition.
The output includes all 18 coordinate vectors and sign tests for each
lattice.  PARI independently also returned `bnfcertify=1`, ordinary class
group `[3]`, and narrow class group `[3]` in the current session.

## Necessary manuscript adjustment

The original sentence that trace-lattice enumeration handles every
class-number-one field except 49,1125,5125,6125 is not supported by the
actual run.  A minimal accurate replacement, using the
separate arithmetic outcomes, is:

> We enumerate the fields within the bounds of
> Proposition~\ref{prop:c11-finite-boundary}, retaining each isomorphism
> class.  Local factorization tests and certified finite Euler products
> exclude all but the fields recorded in the computational supplement.
> Exact enumeration in the narrow-component lattices then excludes fifteen
> further fields, including the field of discriminant 3969.  The remaining
> discriminants are 49, 169, 361, 725, 1125, and 5125; we treat these by the
> dyadic coefficient bounds and the relevant cuspidal-space computations.

This paragraph describes actual computations.  The subsequent treatment
must state the verified dyadic and space results for all six remaining
discriminants.  The computational script proves the finite arithmetic
inputs; it does not formally verify the manuscript's coefficient lemma,
even-class-number argument, or analytic global bounds.

## Resolution in the free package

The three unresolved_search rows record the bounded lattice search alone.
The separate cubic_weight_two_genus and quartic_725_weight_two certificates
now prove S_2=0 for 169,361,725. Their dyadic bounds exclude all other
source weights. Discriminant1125 is closed by S_2=0 and dim S_5>=5 from
quartic_1125_sufficient.sage. Discriminants49 and5125 are treated as
explained in the dyadic note and the manuscript updates. Thus no lattice
residual is an unresolved endpoint of this combined package.
