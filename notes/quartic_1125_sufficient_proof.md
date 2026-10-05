# A free replacement for the D = 1125 endpoint computation

The accompanying executable `quartic_1125_sufficient.sage` establishes the
finite arithmetic inputs below.  It uses only Sage/PassageMath and exact
number-field arithmetic.  It proves the sufficient statements

\[
S_2(1)=0,\qquad \dim S_5(\varepsilon)\geq5>1,
\]

where `epsilon` is the unique conductor-one totally odd character of the
field.  It does **not** claim the old transcript's exact values
`dim S_3=2` and `dim S_5=6`; neither exact value is needed for the product
exclusion once the dyadic bound leaves only source weights 2 and 3.

## Field and central character

Let
\[
F=\mathbb Q(a),\qquad a^4-a^3-4a^2+4a+1=0.
\]
The script verifies `D_F=1125`, `h_F=1`, `h_F^+=2`, and that all three
fundamental units have norm +1.  Thus the unit-signature image is exactly
the even-sign subspace: its rank is 3, and all its vectors have even
parity.  The nontrivial narrow character `epsilon` is therefore totally
odd.  It is the only character compatible with odd parallel weight.

A finite central idele has associated ordinary ideal `(b)` because
`h_F=1`. Since every base-field unit has norm one,
`epsilon((b)) = sign(Norm(b))` is independent of the generator. The
normalized central action is derived explicitly below; this is what
locates the five-dimensional contribution in the required character
sector. The invariant average alone would not identify that sector.

## A uniform stabilizer bound

For a maximal order `O` in any totally definite quaternion algebra over
this quartic field, the reduced unit group `Gamma=O^*/O_F^*` is finite and
embeds into `SO(3)` at a real place.  If a projective element has order
`m`, its trace `t` and norm `n` satisfy

\[
\frac{t^2}{n}-2=\zeta_m+\zeta_m^{-1}
\]

after choosing a primitive eigenvalue ratio.  Consequently the real
cyclotomic field has degree at most 4 and `phi(m)<=8` (the orders 1 and 2
cause no difficulty).  The elementary inequality `phi(m)^2>=m/2` reduces
the check to `m<=128`.  It follows prime by prime from
`phi(p^e)^2/p^e=p^{e-2}(p-1)^2`: every odd-prime factor is at least 1,
and the possible factor for `p=2,e=1` is 1/2.  Exact enumeration gives
`m<=30`.  A cyclic reduced unit group therefore has order at most 30,
a dihedral group at most 60, and the other possibilities have orders
12, 24, and 60.  Hence every stabilizer has order at most 60.

The finite-group classification and finiteness statements used here are
Voight, *Quaternion Algebras*, Proposition 11.5.2 and Lemma 26.5.1:
<https://link.springer.com/chapter/10.1007/978-3-030-56694-4_11> and
<https://link.springer.com/chapter/10.1007/978-3-030-56694-4_26>.

## The quaternion ideal class number is two

Take the quaternion algebra unramified at every finite prime and ramified
at all four real places, and a maximal order in it.  The mass formula gives

\[
M=\frac{2D_F^{3/2}\zeta_F(2)}{(2\pi)^8},\qquad
\#\operatorname{Cls}(O)\leq60M.
\]

The script verifies that 2 is inert and 3 has one prime of residue degree
2 and ramification index 2.  Using these two Euler factors explicitly and
bounding all other rational-prime factors by those of `zeta(2)^4` yields

\[
\zeta_F(2)
 < \zeta(2)^4\frac{27}{85}\frac{256}{405},\qquad
M<\frac{1125\sqrt{1125}}{826200}<\frac1{20}.
\]

The last comparison is checked by squaring positive quantities and
comparing integers.  Thus the ideal class number is less than 3.
Reduced norm maps the right ideal class set onto `Cl^+(F)`: local norms
are surjective at every finite place because the algebra is split there,
and norms of global quaternion elements are totally positive.  Hence
there are at least two classes, and therefore exactly two.

The mass formula reference is Voight, Main Theorem 26.1.5 (formula
26.1.6), with unit reduced discriminant.  This argument uses a bound and
does not require evaluating `zeta_F(-1)` exactly.

The weight-two quaternionic space is the space of functions on these two
classes.  The two distinct characters obtained by composing reduced norm
with the two narrow class characters already span it.  Removing this
noncuspidal subspace leaves `S_2(1)=0` under Jacquet--Langlands.

## A five-dimensional contribution in weight five

Put `tau=a^3-3a+1`, so `tau^2-tau-1=0`.  In
`B=(-1,-1)_F`, with `i^2=j^2=-1` and `ij=-ji`, let

\[
z=\frac{\tau+(\tau-1)i+j}{2},\qquad
O=O_F+O_Fi+O_Fz+O_Fiz.
\]

The script checks all 16 products of the displayed basis, showing that
this is an order, and checks that the reduced-trace Gram determinant is
-1.  Its discriminant is therefore the unit ideal; the order is maximal
and the algebra is unramified at every finite prime.  At the real places
it is Hamilton's algebra, as required.  This is also the scalar extension
of Voight's icosian order in (11.5.8).

Exact breadth-first enumeration shows that `i` and `z` generate precisely
120 distinct elements, each of reduced norm one.  Their image in the
reduced unit group has order 60.  The uniform upper bound above makes
this the full reduced unit group; in particular no larger stabilizer
has been overlooked.

For parallel weight five, restrict the coefficient representation to
these norm-one units.  It is
\[
W=\bigotimes_{\sigma:F\hookrightarrow\mathbb R}
\operatorname{Sym}^3(\mathbb C^2).
\]
For an element of reduced trace `t`, its character is
\[
\operatorname{tr}(g\mid W)
 =N_{F/\mathbb Q}(t^3-2t).
\]
The determinant normalization is trivial on norm-one units, and the
central element -1 acts trivially on the tensor product of four odd
symmetric powers.  The exact 120-element average is
\[
\dim W^{O^*/O_F^*}
 =\frac1{120}\sum_gN_{F/\mathbb Q}(t_g^3-2t_g)
 =\frac{600}{120}=5.
\]
All 120 elements and their character values are saved in the JSON output.
The output also records the distribution `256:2`, `1:88`, `0:30`, the
exact multiplication table, and the full trace Gram matrix.

## The central character of this contribution

Use the algebraic weight action in Dembélé--Voight (7.3) and (7.10).
A scalar `b` acts on the fourfold weight-five coefficient representation
by `Norm(b)^3`. Thus global covariance says

\[
\Phi(bx)=N_{F/\mathbb Q}(b)^{-3}\Phi(x).
\]

Let a finite central idele `c` have associated fractional ideal
`mathfrak b=(b)`. In the inverse-translation convention of (7.14), the
normalized central operator is

\[
(Z_{\mathfrak b}\Phi)(x)
  =(N\mathfrak b)^{-3}\Phi(xc^{-1}).
\]

Full level gives invariance under `c/b` in the finite unit group, so

\[
Z_{\mathfrak b}\Phi
 =\frac{N_{F/\mathbb Q}(b)^3}{|N_{F/\mathbb Q}(b)|^3}\Phi
 =\varepsilon(\mathfrak b)\Phi.
\]

This formula applies to the whole weight-five quaternionic module and
therefore to the principal-class summand. The executable checks both
central sign cosets, using `b=a+1` of norm `-5`, together with rational
scalars and all fundamental units. For example, `b=a+1` has algebraic
scalar action `(-5)^3=-125`, ideal norm `5`, and normalized central action
`-125/5^3=-1`. In weight two the exponent is zero, giving trivial central
character.

In the definite quaternionic description, the principal right ideal
class contributes this invariant space as a direct summand.  The other
class can only add to the dimension.  In weight five no one-dimensional
norm-character representations occur, so this is a cuspidal contribution.
Jacquet--Langlands therefore gives
`dim S_5(epsilon)>=5>1`.

The companion `verify_quaternion_central_character.py` reconstructs the
field and quaternion arithmetic using only rational numbers in Python's
standard library. It independently checks all 120 group elements, all
16 basis products, the trace Gram determinant, the exact character
distribution, and the scalar normalization rows. It reads the newly
generated v2 Sage certificate; the original v1 certificate is rejected
because it lacks the central-action evidence.

This finite verification does not implement the mass formula, the
surjectivity of reduced norm on ideal classes, or Jacquet--Langlands.
Those are the cited mathematical inputs. The ordinary/narrow class
numbers and completeness of fundamental units are obtained by Sage with
proof mode enabled, and their signature vectors are recorded using exact
real-algebraic embeddings.

Run the arithmetic certificate and independent verifier in that order:

    sage -python quartic_1125_sufficient.sage --output rerun/quartic_1125_sufficient.json
    python verify_quaternion_central_character.py --data-dir rerun --output verification_case11/quaternion_central_character.json

For the quaternionic direct-sum description and its relation with
Hilbert forms, use Dembélé--Voight, *Explicit methods for Hilbert modular
forms*, sections 7--8, corrected May 30, 2024:
<https://jvoight.github.io/articles/hmf-crm-bcn-053024.pdf>.

## Suggested replacement for the D = 1125 paragraph

> At discriminant 1125, the dyadic inequality excludes source weights at
> least four.  The finite quaternionic calculation in the supplement gives
> `S_2(1)=0` and `dim S_5(epsilon)>=5`.  Indeed, the unramified definite
> quaternion algebra has two ideal classes, whose norm characters exhaust
> its weight-two module; an explicit maximal icosian order contributes
> five invariant vectors in weight five.  Thus no weight-two source exists,
> and every weight-three source has a target of dimension greater than one.

Only these sufficient conclusions should replace the unavailable Magma
dimension transcript.  The old exact value 6 is not used or certified by
this replacement.
