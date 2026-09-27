#!/usr/bin/env python3
"""Run the free arithmetic certificates, preserving the supplied outputs."""
import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoints-only", action="store_true")
    args = parser.parse_args()
    if sys.flags.optimize:
        raise RuntimeError("Assertions must remain enabled; do not use -O.")
    import sage.all
    out = ROOT / "rerun"
    out.mkdir(exist_ok=True)
    def output(name):
        return ["--output", str(out / name)]
    tasks = [] if args.endpoints_only else [
        ("e3_cusp_check.sage", []),
        ("e2_square_branch_check.sage", []),
        ("e2_high_degree_check.sage", []),
        ("e2_candidate_enumeration.sage", output("e2_candidates.json")),
        ("e2_allclass_balanced.sage", output("d3969_certificate.json")),
        ("e2_indecomposable_search.sage",
         ["--input", str(out / "e2_candidates.json"), *output("e2_lattices.json")]),
        ("e2_dyadic_bound.sage", output("e2_dyadic_output.json")),
        ("rq_certificate.py", output("rq_certificate.json")),
        ("rq_e2_exact.sage", []),
    ]
    tasks += [
        ("cubic_weight_two_genus.sage", output("cubic_weight_two_genus.json")),
        ("quartic_725_weight_two.sage", output("quartic_725_weight_two.json")),
        ("quartic_1125_sufficient.sage", output("quartic_1125_sufficient.json")),
    ]
    scope = "endpoints" if args.endpoints_only else "all_arithmetic"
    status = {"scope": scope, "completed": False, "tasks": [],
              "magma_required": False,
              "theoretical_inputs": "Read notes/; D21 and D12 are textual proofs."}
    status_path = out / (scope + "_status.json")
    for name, options in tasks:
        print("RUN", name, flush=True)
        log = out / (Path(name).stem + ".txt")
        with log.open("w") as handle:
            result = subprocess.run([sys.executable, str(ROOT / name), *options],
                                    cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
        accepted = result.returncode == 0
        if name == "e2_indecomposable_search.sage" and result.returncode == 2:
            data = read(out / "e2_lattices.json")
            accepted = (data["search_complete"] and data["all_required_degrees"]
                        and sorted(r["discriminant"] for r in data["unresolved_fields"])
                        == [169, 361, 725])
        status["tasks"].append({"file": name, "exit_code": result.returncode,
                                "accepted": accepted, "log": log.name})
        status_path.write_text(json.dumps(status, indent=2) + "\n")
        if not accepted:
            print("Not completed: inspect", log)
            return 1
    cubic = read(out / "cubic_weight_two_genus.json")
    assert {r["discriminant"] for r in cubic["rows"]} == {169, 361}
    assert all(r["full_level_parallel_weight_two_cusp_dimension"] == 0
               for r in cubic["rows"])
    quartic = read(out / "quartic_725_weight_two.json")
    assert quartic["full_level_parallel_weight_two_cusp_dimension"] == 0
    large = read(out / "quartic_1125_sufficient.json")
    assert large["conclusions"]["dim_S2_trivial_character"] == 0
    assert large["conclusions"]["dim_S5_totally_odd_character_at_least"] >= 5
    if not args.endpoints_only:
        enumeration = read(out / "e2_candidates.json")
        assert enumeration["covers_all_required_degrees"]
        assert len(enumeration["fields"]) == 771
        lattices = read(out / "e2_lattices.json")
        assert sum(r["status"] == "excluded_by_lattice_certificate"
                   for r in lattices["fields"]) == 15
    status["completed"] = True
    status["new_arithmetic_endpoints_closed"] = [169, 361, 725, 1125]
    status_path.write_text(json.dumps(status, indent=2) + "\n")
    print("FREE_ENDPOINTS_PASSED" if args.endpoints_only
          else "FREE_ARITHMETIC_SUITE_PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
