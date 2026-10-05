# Real quadratic reduction and local coefficient certificates

`rq_certificate.py` uses Python integers and exact fractions. Its output
records every norm and ideal witness. `rq_e2_exact.sage` independently checks
the number-field arithmetic in Sage and verifies the coefficient identities
symbolically. `verify_d21_local.py` checks the complete discriminant-21
coefficient contradiction and the discriminant-40 indecomposability argument.

The reduction covers all 26 positive fundamental discriminants allowed by
the strict ramified and inert inequalities D^3 < 600^2 and D^3 < 2160^2.
Minkowski's bound and explicit principal-ideal relations give the class-number
upper bounds. For D=40,60,85,165, the totally real biquadratic extension with
quadratic discriminants 5, D/5, and D gives the lower bound two: its
discriminant is D^2, so it is unramified everywhere over the quadratic field.
Norm-minus-one unit witnesses or congruence obstructions determine the narrow
class numbers. The generalized Bernoulli formula supplies all 15 rows of the
manuscript's reciprocal-constant table.

After integrality, the ramified lower bound, and the discriminant-40
indecomposability argument, the survivors are exactly 12,21,24,28,69,77.
For D=12,24,28, the ramified coefficient identity and central-character
parity exclude every source weight at least two.

For inert 2 and D>16, all ordered totally positive decompositions of 4 are
(1,3),(2,2),(3,1). Ideal divisor sums and the Hecke recurrences yield

    3*h_(2)+h_(3) = alpha-f_(3)-(15/alpha)*4^(ell-1).

At D=69 and 77, alpha=2 and f_(3)=13 and 10. With t=ell-2>=0,
the required magnitude is strictly greater than 30*4^t, while the
Ramanujan upper bound is at most 12*2^t+9*3^t <=21*3^t.
The checks 30>21 and 4>3 prove the contradiction for every integer t>=0.
The weight-two values -41,-38 and sharper bounds 21,18 are also recorded.

At D=21, the prime above 3 has generator 1+w of norm -3, where
w=(1+sqrt(21))/2. The coefficient formula excludes every even source weight
and every odd source weight at least five by an induction in steps of two.
For the only remaining weight three, the six ordered decompositions of 5
produce coefficient difference -1188. Indecomposability of the prime
exponents above 5 and Hecke multiplicativity instead force zero. No
cusp-space dimension or Hecke determinant is used. The proof is recorded in
`notes/d21_local_free.tex` and its finite arithmetic is checked by
`verify_d21_local.py`.

At D=40, the verifier checks all 18 possible integer summands of
nu=7+2sqrt(10), finds no decomposition, and certifies
(3,sqrt(10)-1)^2=(nu) using exact lattice indices and ideal generators.
The manuscript's Hecke recurrence then supplies the contradiction.

These arithmetic certificates concern source weights at least two. The
manuscript's endpoint pairing and twisting argument separately exclude source
weight one, so no weight-one endpoint is left open in these notes.

Run from the repository root:

    python3 rq_certificate.py --output rerun/rq_certificate.json
    sage -python rq_e2_exact.sage
    python3 verify_d21_local.py --output verification_case11/d21_d40_local_result.json
