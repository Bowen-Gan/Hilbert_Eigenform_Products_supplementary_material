# Finite field screening and component-lattice certificates

This note records the finite arithmetic inputs used in the manuscript's
Eisenstein-weight-two section. The source programs, input tables and their
hashes, machine-readable outputs, and the final run status are supplied
separately. The lattice search's incomplete-success flag is deliberately
preserved: separate dyadic or cusp-space arguments are needed for six
fields, even when the prescribed lattice search has finished.

## Execution

Run from the repository root in a SageMath Python environment:

```sh
sage -python e2_candidate_enumeration.sage --output rerun/e2_candidates.json
sage -python e2_indecomposable_search.sage --input rerun/e2_candidates.json --output rerun/e2_lattices.json
sage -python e2_allclass_balanced.sage --output rerun/d3969_certificate.json
python verify_lattice_witnesses.py --data-dir rerun --output verification_case11/lattice_verification.json
```

The last command uses only Python's standard library. It freshly verifies
supplied certificates; it does not claim to regenerate their Sage outputs.
The Sage lattice search returns exit code 2 when the bounded search leaves
169, 361 and 725 without lattice certificates. This is the recorded search
outcome, not a failed algebraic assertion. `search_complete` describes
completion of the prescribed search, while
`all_nonexception_candidates_certified` describes its narrower outcome.

## Complete field input

After the manuscript's even-ordinary-class-number exclusion, the necessary
condition `n*h <= 14` permits odd `h>1` only in degrees 3 and 4. Thus the
finite input uses the following inclusive, outward-rounded bounds:

| Degree | Discriminant endpoint | Bound used |
|---|---:|---|
| 3 | 4218 | General reciprocal-constant bound |
| 4 | 68133 | General reciprocal-constant bound |
| 5 | 209508 | Class number one bound |
| 6 | 2429064 | Class number one bound |

The endpoints are checked with outward interval arithmetic against
`D^3 < (2*pi^4/3)^(2*n)` or `D^3 < (4*pi^2)^(2*n)`, respectively.
Class-number-one fields in degrees 3 and 4 must also satisfy the sharper
bound. Degrees at least seven are handled by a separate theoretical and
numerical bound, not this enumeration.

The current complete input is taken from the compressed Bordeaux totally
real field tables. Their identities and hashes are recorded. Every field
is retained by its defining polynomial, rather than by discriminant alone;
for example, both quartic fields of discriminant 64525 and both cubic
fields of discriminant 3969 are present.

There are 772 input fields: 143 cubics, 552 quartics, 37 quintics and 40
sextics. Under the recorded order of arithmetic tests the outcomes are:

| Recorded outcome | Number |
|---|---:|
| Ordinary class number even | 15 |
| Sharper class-number-one discriminant bound | 521 |
| `(2)` not inert | 135 |
| Small-prime norm criterion | 43 |
| Reciprocal constant cannot be an integer | 37 |
| Requires a lattice or separate arithmetic certificate | 21 |

Earlier Sage-generated input had omitted the quartic field
`x^4 - 19*x^2 - 24*x + 16` of discriminant 65808. Its ordinary class
number is one, so the sharper bound already excludes it; the complete
input also records its actual dyadic ramification. Adding it changes the
input count and one screening count but leaves the 21 surviving fields
unchanged. The only survivor of odd ordinary class number greater than one
is the cubic field `x^3 - 21*x - 35`, of discriminant 3969 and class number
three.

## Finite Euler filter

For cutoff `B=199`, local factor degrees give an exact rational finite
Euler product `P_B`. The strict bounds

`P_B < zeta_F(2) < P_B*(200/199)^n`

follow by bounding each omitted rational-prime factor by
`(1-p^-2)^(-n)` and using
`product_(m>199) (1-m^-2)^(-1) = 200/199`.
When `h=1`, the functional equation gives

`|alpha| = (4*pi^2)^n / (D^(3/2)*zeta_F(2))`.

Its certified interval is intersected with the integers. An empty
intersection proves nonintegrality. The transcript contains local factor
degrees, the exact product and tail bound, and outward numerical endpoints.
No rational special value is guessed or reconstructed. This filter
excludes, among others, discriminants 1369 and 6125.

## Exact all-component lattice calculation

For the integral ideal `A=P` or `P^2`, let `m` generate `A intersect Z` and
put `t=(m)*A^-1`. Then `t` is integral and `A=(m)*t^-1`. The pair `(m,t)`
can be used whether or not its narrow class is principal. Positive
scaling of both exponent and component lattice preserves decompositions;
there is no search through units and no assumption that all classes are
narrow principal.

The trace Gram matrix is positive definite. LLL is used only to make its
integer basis more convenient, and the unimodularity of the change of
basis is checked exactly. If `b_j^dual` is trace dual to `b_j`, then
`0 << x << m` implies strict bounds on `Tr(x*b_j^dual)`. Outward rational
interval endpoints, rounded away from the permitted range, give a finite
closed integer box containing every possible summand. Each coordinate
vector in the box is tested. Sage interval signs are accepted only when
decisive; otherwise exact algebraic signs are used. A complete empty
list of positive decompositions certifies indecomposability. A bounded
prime search that fails to find such a box records decomposition witnesses
and does not certify exclusion.

The fifteen successful prime certificates are:

| Degree | Discriminant | Prime norm | First box | Second box |
|---|---:|---:|---:|---:|
| 3 | 81 | 3 | 18 | 64 |
| 3 | 257 | 3 | 12 | 48 |
| 3 | 321 | 3 | 18 | 48 |
| 3 | 697 | 5 | 12 | 64 |
| 3 | 1257 | 3 | 18 | 18 |
| 3 | 1489 | 7 | 12 | 80 |
| 3 | 3969 | 3 | 18 | 18 |
| 4 | 2525 | 5 | 54 | 144 |
| 4 | 4205 | 5 | 36 | 256 |
| 4 | 4525 | 5 | 54 | 144 |
| 4 | 8069 | 5 | 54 | 144 |
| 4 | 16317 | 5 | 54 | 72 |
| 5 | 14641 | 11 | 1024 | 1024 |
| 5 | 38569 | 7 | 324 | 576 |
| 6 | 371293 | 13 | 1296 | 3072 |

The thirty boxes contain 8776 coordinate vectors. Machine-readable
certificates record maximal-order and ideal bases, trace Gram matrices,
coordinate bounds and point counts. The independent Python verifier checks
ring and ideal arithmetic with rational coefficients and checks every
embedding sign using Sturm-isolated roots and outward rational intervals.
If an interval cannot determine the sign of a nonzero field element, it
is refined until decisive. Termination follows because the defining
polynomial is irreducible and a nonzero reduced element cannot vanish at
one of its conjugates.

The Python verifier checks that each supplied order discriminant equals
the externally reported field discriminant. Except at discriminant 3969,
this identification of the maximal order relies on the field input; it is
not a separate maximal-order computation. It does not independently
reproduce the full field enumeration, class groups, Fourier coefficient
lemma, or endpoint Rankin--Selberg argument.

## Discriminant 3969

The standalone Sage calculation verifies the field discriminant, ordinary
and narrow class groups `C3`, `P=(3,a+1)`, `Norm(P)=3`, and `(3)=P^3`.
For `P` the component lattice is `P^2`, and for `P^2` it is `P`; the
exponent is 3 in each. Each recorded reduced basis gives the containing
box `[0,1] x [-1,1] x [-1,1]`, with 18 coordinate vectors and no positive
decomposition. All 36 sign records are included. The Python verifier also
uses Dedekind's index criterion at 3 and 7 to establish maximality here
independently.

## Separate remaining arguments

The six remaining discriminants are 49, 169, 361, 725, 1125 and 5125.
The search's `unresolved_search` records at 169, 361 and 725 mean only that
rational primes through 43 give no lattice certificate. These cases are
resolved by the separate dyadic and weight-two-space calculations;
1125 also uses a dimension calculation with the prescribed central
character. Their outputs must be combined with the lattice certificates.
No one of these six fields is relabelled as a lattice exclusion.
