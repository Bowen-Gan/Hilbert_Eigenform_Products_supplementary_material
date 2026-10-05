# Supplementary Material: Eigenform Product Identities for Full-Level Hilbert Modular Forms

This repository provides SageMath programs, independent Python verifiers, complete finite input tables, arithmetic certificates, execution records, and mathematical notes for *Eigenform Product Identities for Full-Level Hilbert Modular Forms*. The sources are synchronized with the manuscript revision checked on 5 October 2026. Magma is not required.

## Complete reproduction

From the repository root, use the Python interpreter belonging to a SageMath installation:

```bash
sage -python run_all.py
```

In a conda environment where `import sage.all` succeeds, the equivalent command is:

```bash
python -u run_all.py
```

The driver executes twelve arithmetic tasks and ten independent verification tasks. Successful completion of all twenty-two tasks in this invocation prints:

```text
FULL_CERTIFICATES_PASSED
```

The result is recorded in `verification/full_run_status.json` and summarized in `STATUS.json`. Inspect the actual process exit code, the success marker, and `completed: true` together. A failed new run invalidates the current complete-run status; older results alone do not establish that the changed sources passed.

The driver records the actual Python/Sage versions, times, commands, exit codes, source and input hashes, and generated artifact hashes. Complete field-table inputs are included in the source snapshot. Assertions must remain enabled: do not use `-O`, `-OO`, or `PYTHONOPTIMIZE`.

The latest generated results are in `rerun/` and `verification_case11/`. Same-name outputs are replaced on rerunning. Redundant earlier `logs/` and top-level generated certificates are omitted from this distribution. Git history retains earlier revisions.

## Arithmetic programs

All eleven existing Sage files are retained. Their roles are listed below.

| Program | Mathematical calculation |
| --- | --- |
| `e3_cusp_check.sage` | Weight-three and weight-four exceptional-field arithmetic, explicit unit norms, character parity, and exact Dirichlet special values. |
| `e2_square_branch_check.sage` | Complete quartic input in the square branch and its exclusions. |
| `e2_high_degree_check.sage` | Septic input and higher-degree local exclusions and Euler bounds. |
| `e2_candidate_enumeration.sage` | Complete degree-three through degree-six field inputs; maximal orders, class numbers, local decompositions, and reciprocal-constant screening. |
| `e2_allclass_balanced.sage` | The all-narrow-component certificate at discriminant 3969. |
| `e2_indecomposable_search.sage` | Trace-dual coordinate boxes and certified indecomposability for the fifteen lattice exclusions. |
| `e2_dyadic_bound.sage` | Exact dyadic decompositions, Euler-product bounds, and integer-weight cutoffs for the six final fields. |
| `rq_e2_exact.sage` | Sage cross-check of the quadratic reduction and coefficient algebra for every relevant integer weight. |
| `cubic_weight_two_genus.sage` | CM and genus arithmetic for discriminants 169 and 361. |
| `quartic_725_weight_two.sage` | CM, mass, and weight-two arithmetic at discriminant 725. |
| `quartic_1125_sufficient.sage` | Weight-two vanishing and the five-dimensional contribution with the specified totally odd central character at discriminant 1125. |

`rq_certificate.py` is the twelfth arithmetic task. It independently performs the quadratic reduction using exact standard-library arithmetic; it is also rerun in the independent suite.

The field reduction now starts from **772 fields**, with degree counts **143, 552, 37, 40**. The formerly omitted quartic polynomial is `x^4 - 19*x^2 - 24*x + 16`, with field discriminant **65808**, polynomial discriminant **4211712**, and power-basis index **8**. Its dyadic ramification is recorded independently. The **21 surviving candidates**, **15 lattice exclusions**, and six final discriminants **49, 169, 361, 725, 1125, 5125** are unchanged.

The full compressed Bordeaux tables and their provenance/hashes are in `data/number_field_tables/`. `number_field_table_inputs.py` reads these inputs without requiring network access. Completeness of the original tables is an external mathematical input; it is not inferred from a Sage enumeration-complete flag. Documented duplicate source rows outside our bounds do not alter the required input set.

## Independent verifiers

These programs use only the Python standard library.

| Program | Independent check |
| --- | --- |
| `verify_field_enumeration.py` | Complete source inputs, all 772 computed field records, integral bases/discriminants, local arithmetic, screening intervals, and the missing field's dyadic ramification. |
| `verify_lattice_witnesses.py` | Rational root isolation, trace-dual boxes, exact lattice arithmetic, all thirty boxes, and adaptive sign refinement. |
| `verify_dyadic_bounds.py` | Local factorization and exact rational bounds for the six dyadic cutoffs, including monotonicity at integer weights. |
| `verify_special_values.py` | Generalized Bernoulli sums and exact special values in the weight-three, weight-four, and dyadic calculations. |
| `verify_d21_local.py` | The coefficient contradiction at discriminant 21 and the local indecomposability/prime-square calculations at discriminant 40. |
| `rq_certificate.py` | The exact quadratic reduction and local exclusions. |
| `verify_rankin_local_factor.py` | The generic polynomial identity, character/twist directions, and the correction factor with character `chi^(-1)` at argument `2+2w`. |
| `verify_quaternion_central_character.py` | Exact quaternion arithmetic, the 120-element group and character average, and the scalar action locating the five-dimensional contribution in the epsilon sector. |
| `verify_d12_weight_one.py` | Low-weight Hilbert-series counts, both narrow components, the Dirichlet values at zero, Eisenstein support, and the normalizing scalar 12. |
| `verify_eisenstein_bounds.py` | Certified rational intervals for the finite numerical comparisons in the Eisenstein--Eisenstein proofs. |

An independent-only run checks the existing arithmetic transcripts rather than generating all Sage data again:

```bash
python run_verification.py --data-dir rerun --output-dir verification_case11
```

Its success marker is `CASE11_ADDITIONAL_ARITHMETIC_PASSED`; this is a compatibility name for the expanded independent suite, not the complete-reproduction marker. Individual verifiers accept `--output PATH`; those requiring certificates also accept `--data-dir PATH`. Their JSON outputs state what they verify and their external dependencies.

The lattice search returns exit code 2 when precisely 169, 361, and 725 remain unresolved by that search. The complete driver accepts this only after checking the complete search record and then runs the endpoint programs that establish the required weight-two vanishing. This accepted intermediate exit is not treated as an unconditional proof of all exclusions.

## Mathematical scope and notes

The calculations verify finite arithmetic and algebra. Their application uses the manuscript's reductions and cited theorems, including the mass/genus formulas and Jacquet--Langlands correspondence. The specified-character contribution at 1125 is a sufficient lower bound, `dim S_5(epsilon) >= 5`; no obsolete exact dimension claim is retained.

The first normalization annotation is addressed computationally by the existing local Rankin verifier and its character algebra. Identifying the global Mok family, its nonzero normalization, ordinary class-group projection, the actual Petersson endpoint, and boundary nonvanishing remain the analytic arguments in the manuscript. Finite calculations do not replace those arguments.

The manuscript now treats the Eisenstein weight-two case with cuspidal weight one as well. Its final twisting argument is theoretical. Products of a weight-one Eisenstein eigenform with a cuspidal eigenform remain unclassified; the discriminant-12 identity is an example, not a classification of that open case.

`notes/` provides the accompanying mathematical explanations. `MANUAL_CHECK_REPORT.md` records the per-file source review and the distinction between source review, independently replayed certificates, and a fresh complete run. `manuscript_changes.md` lists the few text changes required by the computational synchronization.

## Integrity and citation

`SHA256SUMS.txt` lists every distributed file except itself and generated Python caches. The complete driver refreshes it after success. Documentation or other changes require a new manifest; a checksum identifies bytes and does not certify a computation.

Cite a fixed Git commit of the exact source-and-record version used. The repository is [Bowen-Gan/Hilbert_Eigenform_Products_supplementary_material](https://github.com/Bowen-Gan/Hilbert_Eigenform_Products_supplementary_material). The manuscript's previous fixed citation to `8eaf1af` must be replaced after this update with the new commit, rather than continuing to cite the older 771-field input.
