# Synchronization with the checked manuscript

The mathematical classification is unchanged. The following text details require synchronization with the revised computation package.

1. The complete finite degree-three through degree-six input contains 772 fields (143 cubic, 552 quartic, 37 quintic, 40 sextic), not 771. The extra quartic is defined by `x^4 - 19*x^2 - 24*x + 16`, with field discriminant 65808 and power-basis index 8. Its dyadic factorization excludes it independently. The latest manuscript already includes this correction. Twenty-one candidates and the fifteen lattice exclusions remain unchanged.

2. In the lattice-box table, the row for discriminant 8069 must read

   ```latex
   8069 & $x^4-x^3-5x^2+5x+1$ & 5 & 54 & 144\\
   ```

   The previous final entry 108 disagrees with both stored certificates and independent exact re-enumeration. The displayed total 8776 already agrees with the correct entry 144. This correction changes no exclusion.

3. The 1125 certificate now records the central action explicitly, rather than inferring the required central character from the dimension of one stabilizer-invariant space. Its conclusion is the sufficient bound `dim S_5(epsilon) >= 5`. The latest manuscript already supplies this argument.

4. The local Rankin verification now includes the formal character/twist direction and the ordered inducing pair. It does not compute or assert an explicit global normalization for the Mok family. The first annotation is resolved in the manuscript by its own nonzero unit-ideal coefficient normalization and cited analytic arguments.

5. Replace the fixed bibliography commit `8eaf1af` and its full commit URL by the commit of the updated package. A commit should never be guessed or replaced by a moving branch URL in the manuscript citation.

6. Additional independent checks cover the discriminant-12 weight-one example and the finite numerical comparisons in the Eisenstein--Eisenstein sections. They do not change the dimension formulas or the identity `h_5 = 12 E_1(1,epsilon) h_4` already proved in the manuscript.

The former incomplete/obsolete computation notes and redundant earlier execution transcripts have been superseded. Only the latest complete run should accompany this revision.
