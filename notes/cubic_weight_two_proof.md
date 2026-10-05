# Free certificate: weight two vanishes at discriminants 169 and 361

`cubic_weight_two_genus.sage` writes its arithmetic certificate to
`cubic_weight_two_genus.json` (or the requested output path). It uses exact
number-field and finite-character arithmetic, with proof mode enabled.
It does not use Magma, precomputed Hilbert-form dimensions, or absence of
database records.

For the two cubic fields F below, take the quaternion algebra B/F ramified
at two real places and no finite place, and a maximal order O in B. Both
ordinary and narrow class numbers of F are one. Consequently the full
adelic weight-two space is represented by the single Shimura curve
associated with O. The norm-positive and norm-one arithmetic groups agree
projectively because every totally positive unit is a square. The
Jacquet--Langlands correspondence identifies its holomorphic
differentials with the full-level parallel-weight-two Hilbert cusp space.

Shimizu's area formula and Riemann--Hurwitz give

\[
\dim S_2(F)=g=1+\frac{|\zeta_F(-1)|}{4}
              -\frac{e_2}{4}-\frac{e_3}{3}.
\]

Only elliptic orders 2 and 3 occur here. Indeed, if an elliptic cycle has
order q, the field F(zeta_(2q)) is quadratic over F. Thus phi(2q) divides
6. The inequality phi(m)^2 >= m/2 reduces the possibilities to a finite
check, yielding q=2,3,7,9. The last two would force F to be the real cubic
cyclotomic field of discriminant 49 or 81, respectively, and therefore
do not occur for either field here.

For q=2,3 let K_q=F(i),F(sqrt(-3)), respectively. The ring
O_F[zeta_(2q)] is already maximal. The certificate verifies this by
equality of the absolute field discriminant and the tensor-order
discriminant D_F^2*(-4)^3 or D_F^2*(-3)^3. Therefore there is only one
quadratic order contributing for each q. All finite local embedding
factors equal one, since B and O have discriminant and level (1).
The unit norm index is also one: norms from a CM field are totally
positive, and totally positive units of F are squares. The elliptic-cycle
formula consequently reduces to e_q=h(K_q).

The executed exact computation gives:

| D_F | Defining polynomial | zeta_F(-1) | h(F(i)) | h(F(sqrt(-3))) | g |
|---:|---|---:|---:|---:|---:|
| 169 | x^3-x^2-4x-1 | -1/3 | 3 | 1 | 1+1/12-3/4-1/3=0 |
| 361 | x^3-x^2-6x+7 | -1 | 1 | 3 | 1+1/4-1/4-1=0 |

The special values are finite generalized Bernoulli sums, after exact
Gaussian-period identification with the cubic subfields of conductors
13 and 19. The four sextic CM class numbers are computed with
`class_number(proof=True)`. Their defining polynomials and discriminants
are recorded in the JSON transcript.
For each cubic field, the complete 8-element unit-signature image is
checked and recorded using exact real-algebraic embeddings. This supports
the assertion that totally positive units are squares and the norm unit
index in the elliptic-cycle formula.

It follows that S_2(F)=0 for both fields. Combined with the already
executed dyadic cutoffs and the absence of odd weights, this finishes
both additional lattice survivors.

The area, elliptic-cycle, Riemann--Hurwitz and Jacquet--Langlands formulas
are cited mathematical inputs. The executable checks their finite
arithmetic inputs and substitutions; it does not prove those theorems.
Actual execution versions and results are recorded by the repository runner.

## References

* John Voight, *Shimura curves of genus at most two*, Math. Comp. 78
  (2009), 1155--1172: equations (1), (2), Lemma 2.1 and Proposition
  2.3(a). Author's version incorporating corrections:
  https://jvoight.github.io/articles/shimbound-mcom-fixed-errata.pdf
* Lassina Dembele and John Voight, *Explicit methods for Hilbert modular
  forms*, Elliptic curves, Hilbert modular forms and Galois deformations,
  Birkhauser (2013), 135--198: Theorem 3.9 and its following paragraph;
  Section 5 identifies the weight-two dimension with the Shimura-curve
  genus. Author's corrected version:
  https://jvoight.github.io/articles/hmf-crm-bcn-053024.pdf

Run:

    sage -python cubic_weight_two_genus.sage --output cubic_weight_two_genus.json

This certificate makes no assertion for discriminants 725 or 1125.
