# Free certificate: the discriminant-725 weight-two endpoint

The script `quartic_725_weight_two.sage` writes its complete arithmetic
certificate to `quartic_725_weight_two.json` (or the requested output
path). It proves the required vanishing without
Magma and without using a database dimension or a numerically guessed
special value.

Let F be defined by x^4-x^3-3x^2+x+1. The script checks D_F=725 and
h(F)=h^+(F)=1. Let B/F be the quaternion algebra ramified at all four
real places and at no finite place, and O a maximal order. By the
Jacquet--Langlands correspondence, full-level weight-two Hilbert cusp
forms correspond to the orthogonal complement of the constants in
the functions on Cl(O). Consequently it suffices to prove h(O)=1.

We use the Eichler class-number formula in Kirschmer--Voight,
Proposition 5.1 and Lemma 5.3. Since the quaternion discriminant and
the order level are both (1), all finite embedding factors are one.
Every positive unit of F is a square. For any CM order R, dividing
a unit u by a base-field unit of norm-square root makes all its complex
absolute values equal to one. Kronecker's theorem then makes it a root
of unity. Thus there are no extra half-elliptic unit contributions, and
[R^*:O_F^*] equals half the number of roots of unity in R.
The full 16-element unit-signature image is also checked and recorded
using exact real-algebraic embeddings; no approximate sign determination
is used. Class-number and fundamental-unit calculations run in proof mode.

All possible CM roots are found completely: if F(zeta_(2q))/F is
quadratic, phi(2q) divides 8. The bound phi(m)^2 >= m/2 makes q<=64
sufficient. Exact factorization of these cyclotomic polynomials over
F leaves precisely q=2,3,5. For each, the script constructs the
quadratic CM extension and verifies that O_F[zeta_(2q)] is maximal
by equality between its polynomial discriminant ideal and the relative
field discriminant ideal. All conjugates of the root lie in that
quadratic field, so different quadratic factors do not yield further
extension fields. In each case the absolute CM class number is one,
computed with proof mode enabled.

The class-number formula therefore specializes to

\[
h(O)=M+\frac12\left(1-\frac12\right)
       +\frac12\left(1-\frac13\right)
       +\frac12\left(1-\frac15\right)
    =M+\frac{59}{60},
\]

where M is the mass. We only need a coarse rigorous upper bound for M.
The mass formula and the Dedekind functional equation give

\[
M=2\frac{725^{3/2}\zeta_F(2)}{(4\pi^2)^4}.
\]

The elementary inequalities zeta_F(2)<=zeta(2)^4<16 and pi>3 imply

\[
0<M<\frac{2\cdot725\sqrt{725}\cdot16}{36^4}
     <0.371918.
\]

The final comparison uses exact real-algebraic arithmetic. Hence
59/60<h(O)<1.355252. Since h(O) is an integer, it is one. This proves
S_2(F)=0 and finishes the remaining source weight for discriminant 725.
As a consequence, the same formula gives M=1/60 and zeta_F(-1)=2/15;
these exact values are deductions, not numerical reconstructions.

The finite computation supplies arithmetic inputs for the cited mass and
embedding formulas. It does not implement or independently prove those
theorems or the Jacquet--Langlands correspondence. Actual execution
versions and results are recorded by the repository runner.

## Sources and independent agreement

* Markus Kirschmer and John Voight, *Algorithmic enumeration of ideal
  classes for quaternion orders*, SIAM J. Comput. 39 (2010), 1714--1747:
  equation (5.1), Proposition 5.1, Lemma 5.3 and equation (5.3).
  https://jvoight.github.io/articles/quatideal-fixed-errata-111614.pdf
* The same paper's Proposition 8.1 and Table 8.2 also list this maximal
  order as having class number one. That table is an independent check;
  the executed certificate above uses the formula and freshly computed
  CM arithmetic instead of the table as input.
* The later errata state that the class-number-one list was unaffected
  by the table-enumeration bug, and independently verified:
  https://jvoight.github.io/articles/quatideal-errata-102320.pdf
* Lassina Dembele and John Voight, *Explicit methods for Hilbert modular
  forms*, Definition 3.7 and Theorem 3.9:
  https://jvoight.github.io/articles/hmf-crm-bcn-053024.pdf

Run:

    sage -python quartic_725_weight_two.sage --output quartic_725_weight_two.json
