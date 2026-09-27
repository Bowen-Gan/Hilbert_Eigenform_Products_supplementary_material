# Supplementary Material: Eigenform Product Identities for Full-Level Hilbert Modular Forms

This repository contains executable arithmetic checks, recorded outputs, and mathematical notes supporting the finite computations in *Eigenform Product Identities for Full-Level Hilbert Modular Forms*.

The scripts verify the arithmetic inputs to the arguments identified below. Their interpretation uses the reductions and coefficient identities in the manuscript, together with the theorems cited in the accompanying notes. The package uses SageMath and Python; Magma is not required.

## Recorded verification results

The repository includes the following execution records:

| Record | Reported result |
| --- | --- |
| [`rerun/all_arithmetic_status.json`](rerun/all_arithmetic_status.json) | The complete arithmetic suite finished with `completed: true`; all 12 tasks have `accepted: true`. |
| [`rerun/endpoints_status.json`](rerun/endpoints_status.json) | The three endpoint programs finished with `completed: true`; each has exit code zero. |
| [`rerun/`](rerun/) | Individual program logs and structured arithmetic outputs. |
| [`STATUS.json`](STATUS.json) | A package-level summary, dated 27 September 2026, including the scope of the replacement arguments. |

The endpoint JSON files in `rerun/`, as well as the enumeration, lattice, and dyadic outputs there, record Sage version **10.9**. The corresponding original endpoint JSON files at the repository root record version **10.8.12**. These are separate sets of recorded outputs; the version in `STATUS.json` describes the original set.

The records supplied here document successful completion of the implemented checks. They do not contain a complete record of the execution environment: a run-specific timestamp, operating-system information, Python version, and source revision captured at execution time are not recorded by the current driver. Such information should accompany any newly archived local run; a procedure is given below.

A successful execution of this suite verifies its implemented arithmetic checks. It is not a formal verification of the manuscript or of the cited mathematical theorems.

### The lattice-search exit status

The full-suite record contains one deliberately accepted nonzero exit code. The program `e2_indecomposable_search.sage` returns code `2` when its bounded search leaves unresolved fields. The driver accepts this outcome only after checking that:

- the search is complete over its prescribed input and covers all required degrees; and
- the unresolved fields have exactly the discriminants `169, 361, 725`.

The endpoint programs subsequently verify the vanishing of the weight-two cusp spaces for those fields. Consequently, `unresolved_search` and `Nonexception lattice certification complete: False` in the individual lattice output describe that search alone. The full-suite success condition also requires the separate endpoint checks. All other tasks in the committed full-suite record have exit code zero.

Some endpoint text logs are empty because those programs write their results to JSON files. Their results and acceptance status must be read together with the corresponding JSON output and the driver status file.

## Reproducing the computations

Use a Python environment that provides SageMath. The recorded rerun outputs use Sage 10.9. Run the following commands from the repository root after activating that environment:

```bash
python -c "import sage.all; from sage.version import version; print(version)"
python run_all.py
```

The final success marker for the complete suite is:

```text
FREE_ARITHMETIC_SUITE_PASSED
```

To run only the three endpoint programs:

```bash
python run_all.py --endpoints-only
```

The corresponding success marker is:

```text
FREE_ENDPOINTS_PASSED
```

An endpoint-only run does not execute the other nine programs. Where the Sage launcher supports `sage -python`, that command may be used in place of `python`.

Python assertions must remain enabled: do not use `-O`, `-OO`, or an optimization setting in `PYTHONOPTIMIZE`. The driver rejects an optimized Python invocation.

The driver writes individual logs, JSON results, and a scope-specific status file to `rerun/`. **Files of the same name in that directory are overwritten.** Use a fresh checkout or preserve the supplied `rerun/` directory before reproducing the computations. A status file from an earlier run is not evidence that a later invocation succeeded; check the current process result and its newly produced records.

## Contents and mathematical scope

Here `D` denotes the discriminant of the totally real field. Cusp-space notation and characters follow the manuscript and the linked notes.

| Program or document | Purpose |
| --- | --- |
| [`e3_cusp_check.sage`](e3_cusp_check.sage) | Finite number-field, unit-signature, primitive-character, and special-value checks for Eisenstein weight three. |
| [`e2_square_branch_check.sage`](e2_square_branch_check.sage) | Quartic enumeration, ideal factorizations, and inequalities for the dyadic square branch. |
| [`e2_high_degree_check.sage`](e2_high_degree_check.sage) | Septic enumeration, septic and octic local checks, and finite Euler-product bounds. Completeness of the octic field list uses the table cited in the manuscript. |
| [`e2_candidate_enumeration.sage`](e2_candidate_enumeration.sage) | Enumeration within the degree-dependent bounds supplied by the theoretical reduction. The recorded output contains 771 fields and 21 candidates requiring further lattice or exceptional-case arguments. |
| [`e2_allclass_balanced.sage`](e2_allclass_balanced.sage) and [`e2_indecomposable_search.sage`](e2_indecomposable_search.sage) | Lattice certificates in narrow-class components, including the separate certificate for `D = 3969`. The driver checks that 15 candidate fields receive lattice certificates. |
| [`e2_dyadic_bound.sage`](e2_dyadic_bound.sage) | Arithmetic inputs to the dyadic inequalities and remaining weight restrictions. |
| [`rq_certificate.py`](rq_certificate.py) and [`rq_e2_exact.sage`](rq_e2_exact.sage) | The six-field real-quadratic reduction and local arithmetic checks. The first program uses the Python standard library; the second checks the arithmetic using Sage. |
| [`cubic_weight_two_genus.sage`](cubic_weight_two_genus.sage) | CM class-number and genus computations supporting `S_2 = 0` for `D = 169, 361`; see [the cubic proof](notes/cubic_weight_two_proof.md). |
| [`quartic_725_weight_two.sage`](quartic_725_weight_two.sage) | Quaternion class-number inputs supporting `S_2 = 0` for `D = 725`; see [the quartic proof](notes/quartic_725_weight_two_proof.md). |
| [`quartic_1125_sufficient.sage`](quartic_1125_sufficient.sage) | Arithmetic supporting `S_2(1) = 0` and `dim S_5(epsilon) >= 5` for `D = 1125`; see [the proof and character specification](notes/quartic_1125_sufficient_proof.md). |
| [`notes/d21_local_free.tex`](notes/d21_local_free.tex) | A coefficient argument for `D = 21`, replacing the earlier Hecke-determinant calculation. |
| [`notes/d12_ring_dimensions.md`](notes/d12_ring_dimensions.md) | A derivation from published ring structures for `D = 12`. |

The scripts use exact integer, rational, and number-field arithmetic, with real interval or ball arithmetic where analytic bounds are needed. A displayed decimal approximation is not a substitute for the comparison implemented in the script.

The [`notes/`](notes/) directory records the mathematical interpretation and external theorem dependencies. The [manuscript revision checklist](manuscript_changes.md) identifies passages that must agree with these certificates. In particular, the current package does not certify the former `D = 21` Hecke determinants or the former exact values `dim S_3 = 2` and `dim S_5 = 6` at `D = 1125`. The replacement arguments use the sufficient conclusions specified above.

## Archiving a local run

For a local execution record, retain the command transcript, the newly generated `rerun/` directory, and the environment and source information. The following example uses Bash on Linux or WSL, from an activated Sage environment and a Git checkout of this repository. It creates a separate archive directory for each invocation:

```bash
(
set -e -o pipefail
mkdir -p verification
run_dir=$(mktemp -d "verification/run-$(date -u +%Y%m%dT%H%M%SZ)-XXXXXX")

{
    date --iso-8601=seconds
    uname -srm
    python -c "import sys; from sage.version import version; print('Sage:', version); print('Python:', sys.version)"
    git rev-parse HEAD
    git status --short
    git diff HEAD -- run_all.py '*.sage' rq_certificate.py
    sha256sum run_all.py *.sage rq_certificate.py
} > "$run_dir/environment.txt" 2>&1

if python -u run_all.py 2>&1 | tee "$run_dir/local-run.txt"; then
    run_status=${PIPESTATUS[0]}
else
    run_status=${PIPESTATUS[0]}
fi
printf '\nDriver exit code: %s\n' "$run_status" \
    | tee -a "$run_dir/local-run.txt"
date --iso-8601=seconds > "$run_dir/finished-at.txt"
cp -a rerun "$run_dir/rerun"
printf 'Execution record saved in %s\n' "$run_dir"
exit "$run_status"
)
```

The driver must terminate with exit code zero and print the full-suite success marker. Its newly generated `all_arithmetic_status.json` must report `completed: true`. After saving the record, the subshell above returns the driver's exit code; a failure while creating or copying the archive also causes a nonzero exit. The copied `rerun/` directory is preserved independently of subsequent invocations. If the run fails, the archived directory may contain incomplete results or files left by an earlier invocation; it must not be described as a successful verification.

The source diff is taken against `HEAD`, so it includes both staged and unstaged changes to the tracked program files. For a published verification record, use committed source files and identify the exact commit that was executed. These records supplement the source code and mathematical arguments; they do not replace them.

To archive an endpoint-only run, replace `python -u run_all.py` in the command above with `python -u run_all.py --endpoints-only`. Require the endpoint success marker and a newly generated `endpoints_status.json` with `completed: true`, and describe the record explicitly as an endpoint run. Other files copied from `rerun/` do not establish that the full suite was executed in that invocation.

The files under `verification/` in this example are to be generated on the machine performing the run; they are not included among the existing records.

## File integrity

[`SHA256SUMS.txt`](SHA256SUMS.txt) lists checksums for the distributed source, documentation, and recorded output files, including `rerun/`. It excludes itself and generated Python cache files. Before rerunning or editing the files, the manifest can be checked from the repository root with:

```bash
sha256sum -c SHA256SUMS.txt
```

Checksums identify file contents; they do not by themselves establish that a computation was executed. The run records and mathematical arguments serve different purposes.

When distributing a revised package, update the checksum manifest to match the files being distributed, including any revised documentation or recorded outputs. Changes to this README also require its checksum to be updated.

## Citation

Cite the version of this supplement that accompanies the manuscript under review. A citation should identify the supplement, the repository, and either the full commit identifier or the DOI of a specific archived release. A link to the repository's default branch alone does not identify a fixed version.

The repository is available at:

<https://github.com/Bowen-Gan/Hilbert_Eigenform_Products_supplementary_material>

For a manuscript using `thebibliography`, the following template gives a title-based citation. Replace `FULL_COMMIT_SHA` with the complete identifier of the version used, and `SHORT_COMMIT_SHA` with its abbreviated identifier. Obtain the full identifier from that checkout using `git rev-parse HEAD`. Add the authors of the supplement before the title once its authorship has been specified.

```latex
\bibitem{EigenformCertificates}
\emph{Supplementary Material: Eigenform Product Identities for
Full-Level Hilbert Modular Forms},
SageMath source code, arithmetic certificates, and recorded outputs,
GitHub repository, 2026,
commit \texttt{SHORT_COMMIT_SHA},
\url{https://github.com/Bowen-Gan/Hilbert_Eigenform_Products_supplementary_material/tree/FULL_COMMIT_SHA}.
```

The supplement may be cited in the manuscript as follows:

```latex
The SageMath source code, arithmetic certificates, recorded outputs,
and instructions for reproducing the computations are available
in the supplementary material~\cite{EigenformCertificates}.
```

If a release is archived in Zenodo, cite its specific-version DOI to identify the files used for the manuscript. Include the actual release version and authorship recorded in the archive. A concept DOI refers to the collection of versions and does not identify the particular version used. A `CITATION.cff` file may be added at the repository root to provide consistent citation metadata through GitHub's citation interface.
