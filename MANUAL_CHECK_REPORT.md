# Manual source and certificate review

This report records the source review against the checked 5 October 2026 manuscript. Source inspection, independent certificate replay, and regeneration by Sage are different checks. The authoritative complete-run results are the actual files linked from `verification/full_run_status.json`; they are not inferred from this prose.

## Review sequence

1. All 86 files of the previous GitHub revision `8eaf1af86dfc3ea4a9e1fd32114dfa59054ae105` were matched byte-for-byte to their Git blob identities before editing.
2. Every computation was matched to its mathematical reduction, parity/character conditions, finite input bound, and result used in the manuscript.
3. Each program's arithmetic, field/ideal normalization, rounding direction, completeness argument, error handling, and output interface was read and checked.
4. Complete external field inputs were reconciled with the computed records. The omitted quartic was recomputed independently, and the source-table duplicate rows outside our cutoffs were documented explicitly.
5. Exact standard-library verifiers were used to replay the finite certificates. Actual producer runs were then performed with genuine SageMath, followed by a final complete run from frozen source files.
6. Generated records were checked against their current source/input hashes; JSON schemas, expected counts, all referenced artifacts and logs, and the distribution checksum manifest were checked before publication.

## Every executable source file

| File | Manual checks and resulting behavior |
| --- | --- |
| `cubic_weight_two_genus.sage` | Relative CM discriminants, signatures/class numbers, roots of unity, embedding corrections and genus formula for 169/361. External genus and Jacquet--Langlands theorems are stated. |
| `e2_allclass_balanced.sage` | Integral component lattices, unimodular basis changes, trace-dual coordinate bounds and exact positivity fallback; both 3969 exponents. |
| `e2_candidate_enumeration.sage` | Every complete table row, endpoint ceilings, maximal orders/indices, certified class numbers, local factors and reciprocal-constant screening; 772 inputs and 21 survivors. |
| `e2_dyadic_bound.sage` | Character special values, exact Euler-product/tail directions, outward square-root bounds, complete inverse-Gram short-vector boxes and strict positivity; integer-weight monotonicity. |
| `e2_high_degree_check.sage` | Complete septic input, explicit externally sourced octic list, actual maximal-order dyadic factors, small-prime exclusions and finite Euler products. |
| `e2_indecomposable_search.sage` | All 21 candidate/status records; component representatives and exhaustive prime/prime-square boxes for the fifteen successful certificates; failed bounded attempts remain unresolved. The driver accepts exit 2 only for exactly 169/361/725. |
| `e2_square_branch_check.sage` | All 17 bounded quartics, actual ramification/residue degrees, proof-enabled class numbers and exact square-branch endpoint. |
| `e3_cusp_check.sage` | Full small-field input, explicit negative-norm units, character parity/primitive conductors, weight-three/four Bernoulli products and rational cutoffs. |
| `number_field_table_inputs.py` | Safe literal-list parsing rather than GP execution, source and compressed hashes, coefficient order, cutoff strictness, repeated discriminants and documented source duplicates. |
| `quartic_1125_sufficient.sage` | CM/mass input, full unit-sign image, order closure and trace lattice, 120-element quaternion group, character average and explicit epsilon scalar action. |
| `quartic_725_weight_two.sage` | CM quartic arithmetic, embedding corrections, signatures and weight-two mass/genus calculation. |
| `rq_certificate.py` | Integer Kronecker/Bernoulli arithmetic, Minkowski ideal representatives and principal witnesses, norm-minus-one certificates, all-weight inequalities at 69/77, quadratic survivor list. |
| `rq_e2_exact.sage` | Independent number-field/class-group/ideal verification, exact divisor factorizations, and polynomial identities for every weight rather than finite weight samples. |
| `run_all.py` | All arithmetic and verifier tasks, exact accepted intermediate exit, fresh output generation, source/table snapshots, linking both suites to one invocation, failure invalidation and current summary. |
| `run_verification.py` | Ten independent checks, input/artifact hashes, exact cross-suite consistency, linkage to the active parent run and rejection of stale success as a complete reproduction. |
| `verify_d12_weight_one.py` | Full versus unit-invariant Hilbert-series counts, principal mixed versus nonprincipal upper component, both L(0) values, support parity and coefficient 12. |
| `verify_d21_local.py` | Complete decompositions, even/odd all-weight coefficient induction, discrepancy -1188, and exact discriminant-40 ideal and indecomposability arithmetic. |
| `verify_dyadic_bounds.py` | Rational finite-field/local and trace computations, all short-vector coordinates, positivity, divisor sums, Euler bounds and six cutoff comparisons. |
| `verify_eisenstein_bounds.py` | Thirty-six strict finite comparisons using exact Machin/Taylor rational intervals; infinite-range induction and monotonicity dependencies remain explicit. |
| `verify_field_enumeration.py` | All 772 table/computed correspondences, rational order/ideal arithmetic, Euler intervals and statuses, and ramification independently detected by Frobenius on the maximal order modulo two. |
| `verify_lattice_witnesses.py` | Exact rational root isolation, adaptive refinement for ambiguous nonzero signs, Gram/basis matrices, all boxes and all decomposition decisions. |
| `verify_quaternion_central_character.py` | Independent rational reconstruction of field/quaternion multiplication, group closure/trace distribution, character average and both central-sign cosets. |
| `verify_rankin_local_factor.py` | Generic polynomial cross multiplication, repeated-root specialization, arithmetic/unitary exponents, chi-inverse correction and ordered twist direction. |
| `verify_special_values.py` | Primitive characters and exact generalized Bernoulli sums for every exceptional special-value table; independent rational comparison of the numerical cutoffs. |

Programs whose certificates use assertions reject optimized Python. The complete driver also rejects optimized execution. Sage API calls unavailable in older Sage were replaced by portable exact arithmetic where necessary.

## Supporting mathematical files

All eight inherited notes were reviewed and synchronized; `notes/number_field_tables.md` was added.

| File | Checked scope |
| --- | --- |
| `notes/cubic_weight_two_proof.md` | Genus/mass interpretation and external correspondence inputs. |
| `notes/d12_ring_dimensions.md` | Aoki/Dursthoff series, component identification, unit invariants, epsilon support and normalization. |
| `notes/d21_local_free.tex` | Full coefficient proof, correct verifier reference, no obsolete determinant input. |
| `notes/dyadic.md` | Current six-field table and exact integer-weight cutoff argument. |
| `notes/lattices.md` | Complete input counts, all-narrow-component method, finite coordinate bounds and exact certificate interpretation. |
| `notes/number_field_tables.md` | Bordeaux source identity, range/counts, duplicate rows and externally supplied completeness. |
| `notes/quadratic.md` | Current all-weight local exclusions and separation from the theoretical weight-one twisting proof. |
| `notes/quartic_1125_sufficient_proof.md` | S2 vanishing and the sufficient S5(epsilon) lower bound with explicit central action. |
| `notes/quartic_725_weight_two_proof.md` | CM, mass and genus arithmetic at 725. |

The README, synchronization notes, status records, raw table provenance and checksum manifest were separately checked for consistency with the actual source/data version. Redundant historical generated files are removed; newly produced JSON and logs are checked individually by their schemas, source/input links and hashes.

## Confirmed arithmetic and transcript changes

- Complete input counts are 143, 552, 37, 40, totaling 772. The additional quartic has field discriminant 65808, polynomial discriminant 4211712 and index 8. Actual dyadic factors are two primes, each with ramification index 2, residue degree 1 and norm 2.
- The 21 survivors, 15 lattice exclusions, six final fields, 30 box sizes and 8776 coordinate vectors are unchanged after regeneration. Some basis matrices have harmless unimodular representation changes between Sage versions.
- The manuscript's D8069 second box is 144, not 108. Both previous transcripts, a fresh Sage run and independent rational enumeration agree; the manuscript's total 8776 already used 144.
- Discriminant 3969 has two boxes of 18 vectors, with all 36 sign rows checked.
- Discriminant 1125 has group order 120 and character distribution 256 (twice), 1 (88 times), 0 (30 times). The average is 5. The scalar with negative norm is explicitly checked: Norm(a+1)=-5, and its normalized central action is epsilon((a+1))=-1.
- The weight-one example has full cusp dimensions S4(1)=2 and S5(epsilon)=1. The principal mixed component has unit-invariant modular dimension 2, not the pre-invariant SL2 dimension 3. Its Eisenstein constant is 1/12; the other component is annihilated.
- The first normalization annotation has a finite algebra check, not a fictitious numerical certificate of an analytic identification. The global kernel and actual endpoint integral remain the manuscript's analytic arguments.

## Additional rejection checks

Exact verifier checks were exercised on corrupted boxes and quaternion records. They reject an incorrect box count, an incorrect negative central action, an altered character value, a duplicated group element, an altered trace Gram entry and the obsolete quaternion schema. Adaptive root refinement was exercised on an initially ambiguous nonzero algebraic sign. These checks test failure behavior as well as successful fixtures.

## Limits of finite verification

The original number-field tables' completeness, certified Sage class-group/maximal-order algorithms, automorphic correspondences, ring descriptions and analytic Rankin--Selberg results are specified inputs. This package does not claim to formalize those theorems. A finite comparison is not presented as a proof of an infinite monotonicity assertion. The currently open Eisenstein-weight-one product classification remains open.

## Earlier checked complete execution

The frozen sources were executed in one invocation, run ID `2e43e364-ee73-4db6-aa1b-d791a6be1e24`, from 2026-10-05T16:07:25.088841+00:00 to 2026-10-05T16:09:16.179525+00:00, with genuine SageMath 9.5. All twelve arithmetic tasks and ten independent checks met their stated acceptance conditions. The actual driver exited with code 0 and emitted `FULL_CERTIFICATES_PASSED`. The intermediate lattice exit code 2 has the explicitly documented unresolved-field meaning.

Before publication, every recorded source, input, output and log hash was rechecked against its file; both suites were linked to this single run. All distributed JSON files were parsed, every program was syntax checked, and two fresh runner failure checks confirmed that an old success cannot survive missing Sage or missing verifier inputs. Their record is `verification/runner_failure_checks.json`.

## Latest distributed execution

The latest local complete execution is recorded in
`verification/full_run_status.json` and summarized in `STATUS.json`.
These records give its run ID, times, software versions and task results.
Its console output is preserved in `verification/local-run.txt`.
