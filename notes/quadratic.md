# Case 11: completed real-quadratic certificates

`rq_certificate.py` uses Python's standard library and exact integers and
fractions. Its JSON output records every witness, not just the final table.
`rq_e2_exact.sage` independently checks the number-field arithmetic in Sage
and verifies the local coefficient identities symbolically.

The Python certificate covers all 26 positive fundamental discriminants
under the paper's ramified/inert bounds. It verifies every class-group upper
bound using prime ideals below the Minkowski bound and explicit principal
ideal relations. For class number two it also records the three fundamental
discriminants of the totally real, everywhere unramified biquadratic
extension, which supplies the lower bound. Unit signatures are certified by
norm-minus-one elements or a congruence obstruction. The complete table of
reciprocal constants is computed by the Bernoulli formula; all 15 rows agree
with the manuscript. The survivors are 12, 21, 24, 28, 69, and 77.

For the inert local calculation, D > 16 proves that all totally positive
decompositions of 4 in the integral basis are (1,3), (2,2), (3,1). The script
computes the ideal-divisor sums at 2 and 3, so the convolution and Hecke
recurrences give

    3 h_(2) + h_(3) = alpha - f_(3) - (15/alpha) 4^(ell-1).

At D=69 and D=77, alpha=2 and f_(3)=13 and 10. At ell=2 this produces -41
and -38, with Ramanujan bounds 21 and 18. These constants are now derived
inside the certificate, rather than supplied in a comment. The dyadic
terminal bounds and the D=21 norm sign are also checked exactly.

The D=40 indecomposability proof remains in the manuscript, where it is
already explicit; the script checks only the ideal and norm witness needed
to connect that proof with the reduction. It does not duplicate that proof.

Run:

    python3 rq_certificate.py --output logs/rq_certificate.json
    sage -python rq_e2_exact.sage > logs/rq_e2_exact.txt

D=21 is now excluded entirely by the coefficient proof in d21_local_free.tex.
No cusp-space dimension or Hecke determinant is needed.
