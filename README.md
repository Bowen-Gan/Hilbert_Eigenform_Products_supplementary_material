# Supplementary Material: Eigenform Product Identities for Full-Level Hilbert Modular Forms

This repository contains SageMath and Python programs, arithmetic certificates, recorded outputs, and mathematical notes supporting the finite computations in *Eigenform Product Identities for Full-Level Hilbert Modular Forms*. Magma is not required.

## Running the complete verification

Use a Python environment in which SageMath is installed and available to Python. Run the following command from the repository root:

```bash
python run_all.py
```

The driver runs the original twelve arithmetic tasks, followed by six additional exact verification tasks. A successful complete run exits with code zero and prints:

```text
FULL_CERTIFICATES_PASSED
```

Use this marker and the corresponding `completed: true` record in `verification/full_run_status.json` as the criterion for a successful complete reproduction. That record is created by the complete driver and links both suites from the same invocation. For subsequent runs, use the latest successful complete run and retain its generated records together.

Outputs and logs are written to [`rerun/`](rerun/) and [`verification_case11/`](verification_case11/). The driver records task outcomes, UTC times, environment information, source-file hashes, and relevant artifact hashes. Existing files with the same names are overwritten, so preserve any outputs that you wish to keep before rerunning. A failed run does not certify completion; older output files alone are insufficient.

Python assertions must remain enabled. Do not use optimized Python execution (`-O`, `-OO`, or `PYTHONOPTIMIZE`).

## Supplied execution records

The repository includes a successful complete run of all eighteen tasks in one invocation of `python run_all.py`, performed on 30 September 2026 using SageMath 10.9 and Python 3.12.14 under WSL2. The run began at 11:59:20 UTC and finished at 12:01:03 UTC. All eighteen tasks were accepted, and the driver reported `FULL_CERTIFICATES_PASSED`.

The complete record is [`verification/full_run_status.json`](verification/full_run_status.json). It links the original twelve-task record, [`rerun/all_arithmetic_status.json`](rerun/all_arithmetic_status.json), and the additional six-task record, [`verification_case11/status.json`](verification_case11/status.json), from that invocation. [`STATUS.json`](STATUS.json) summarizes this latest successful complete run. The corresponding logs and results are included in the two output directories.

The lattice search has one accepted exit code of `2`: its prescribed search leaves precisely the discriminants `169, 361, 725` unresolved. The driver checks this condition, and the subsequent endpoint programs establish the required weight-two cusp-space vanishing. Thus the lattice search's incomplete-exclusion flag must be read together with the endpoint results. All other tasks have exit code zero in the supplied records.

## Computational scope

The original twelve-task suite checks the finite number-field reduction, local arithmetic, lattice exclusions, and remaining cusp-space endpoints:

| Component | Computation |
| --- | --- |
| Eisenstein weight three | Field, unit-signature, character, and special-value checks. |
| Square and higher-degree branches | Quartic and septic enumeration, dyadic ideal factorizations, local inequalities, and finite Euler-product bounds. |
| Candidate enumeration | 771 fields within the degree-dependent bounds, yielding 21 candidates requiring further arguments. |
| Narrow-class lattices | Fifteen candidate-field exclusions, including a separate certificate for discriminant `3969`. |
| Dyadic and quadratic branches | Dyadic weight restrictions, real-quadratic reduction, and local coefficient exclusions. |
| Cusp-space endpoints | Weight-two vanishing for discriminants `169, 361, 725, 1125`, and the sufficient bound `dim S_5(epsilon) >= 5` for `1125`. |

The six additional tasks independently recheck selected arithmetic certificates:

| Program | Exact verification |
| --- | --- |
| [`verify_lattice_witnesses.py`](verify_lattice_witnesses.py) | Fifteen lattice exclusions, thirty exponent boxes, and 8,776 lattice points. |
| [`verify_dyadic_bounds.py`](verify_dyadic_bounds.py) | Six dyadic cutoffs and integer-weight monotonicity. |
| [`verify_special_values.py`](verify_special_values.py) | Four special-value constants; the suite also compares them with the dyadic inputs. |
| [`verify_d21_local.py`](verify_d21_local.py) | The discriminant `21` coefficient contradiction for `ell >= 2`, and discriminant `40` indecomposable-element and prime-square arithmetic. |
| [`rq_certificate.py`](rq_certificate.py) | The quadratic reduction and local exclusions. |
| [`verify_rankin_local_factor.py`](verify_rankin_local_factor.py) | A spherical local polynomial identity and its correction parameter `2 + 2w`. |

The programs use exact integer, rational, and number-field arithmetic, with interval or ball arithmetic for the specified analytic comparisons. Printed decimal approximations do not replace these comparisons.

## Mathematical interpretation

These computations verify the implemented arithmetic. Their application to the classification uses the manuscript's reductions, coefficient identities, and cited theorems. The weight-two Eisenstein–cuspidal checks concern parallel cuspidal weight `ell >= 2`; they do not establish the case `ell = 1`.

The octic list's completeness remains an external table input. The additional lattice verifier does not repeat the full field enumeration or independently identify every maximal order; its output states these dependencies. The local Rankin check does not establish global projection, analytic continuation, or boundary nonvanishing.

For discriminant `1125`, the certified dimension conclusion is the sufficient lower bound above, rather than the earlier exact dimension claims. The six quadratic survivors `12, 21, 24, 28, 69, 77` are intermediate candidates. Earlier discriminant `21` Hecke determinants are not certified by this repository.

The [`notes/`](notes/) directory explains supporting arguments and external dependencies. Some notes and [`manuscript_changes.md`](manuscript_changes.md) retain earlier proof routes and have not yet been synchronized with the current Case 11 proof.

## Integrity and citation

[`SHA256SUMS.txt`](SHA256SUMS.txt) lists hashes of the distributed files, excluding itself and generated Python caches. The complete driver refreshes this manifest after success. Documentation or artifact changes also require a manifest update. Checksums identify file contents; they are not execution certificates.

Cite the supplement using the full Git commit of the version used, or the specific-version DOI of an archived release. The repository is available at [Bowen-Gan/Hilbert_Eigenform_Products_supplementary_material](https://github.com/Bowen-Gan/Hilbert_Eigenform_Products_supplementary_material). A fixed revision identifies the corresponding source and records.
