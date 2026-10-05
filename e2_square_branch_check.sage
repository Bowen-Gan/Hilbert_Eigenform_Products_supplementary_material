"""Exact quartic certificate for the square-ramified dyadic branch.

Read every complete Bordeaux quartic-table row through 4446 (an inclusive
round-up of the strict analytic bound 4445.506), select exactly (2)=p^2
with N(p)=4, and recompute ordinary class numbers with proof enabled.
Table completeness and the modular-form argument are external inputs.
Run: sage -python e2_square_branch_check.sage --output rerun/e2_square_branch.json
"""
from sage.all import *
from sage.version import version as SAGE_VERSION
from number_field_table_inputs import table_rows, load_manifest
import argparse
import json
from pathlib import Path

proof.all(True)
if not __debug__:
    raise RuntimeError("verification assertions require Python without -O/-OO")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="rerun/e2_square_branch.json")
    args = parser.parse_args()
    R = RealIntervalField(256)
    cutoff = ((2*R.pi()**4/3)**4/60)**(QQ(2)/3)
    assert cutoff < R("4445.506") < 4446
    fields = table_rows(4, 4446)
    assert len(fields) == 17
    square_rows = []
    checked_rows = []
    for source in fields:
        polynomial = PolynomialRing(QQ, "x")(source["polynomial"])
        F = NumberField(polynomial, "a")
        assert F.degree() == 4 and F.signature() == (4, 0)
        assert F.discriminant() == source["discriminant"]
        factorization = list(F.ideal(2).factor())
        local_data = [{"norm": int(P.norm()), "ramification_index": int(e),
                       "residue_degree": int(P.residue_class_degree())}
                      for P, e in factorization]
        row = {"discriminant": source["discriminant"],
               "polynomial": source["polynomial"],
               "table_line": source["table_line"],
               "dyadic_factorization": local_data,
               "square_ramified": len(factorization) == 1
                   and factorization[0][1] == 2
                   and factorization[0][0].norm() == 4}
        if row["square_ramified"]:
            class_number = ZZ(F.class_group(proof=True).order())
            assert class_number == source["table_ordinary_class_number"] == 1
            row["ordinary_class_number"] = int(class_number)
            row["class_number_proof"] = True
            square_rows.append(row)
        checked_rows.append(row)
    expected = [1600, 2000, 2624, 3600, 4400]
    assert [row["discriminant"] for row in square_rows] == expected
    assert R(4000)/R.pi()**8 > QQ(4)/15
    payload = {"schema": "case11-square-ramified-quartic-v1",
               "sage_version": SAGE_VERSION,
               "input_provenance": load_manifest()["tables"]["4"],
               "inclusive_enumeration_cutoff": 4446,
               "strict_analytic_cutoff_interval": str(cutoff),
               "complete_table_input_count": len(fields),
               "fields": checked_rows,
               "square_ramified_discriminants": expected,
               "rational_comparison_verified": True}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    print("quartic fields checked from complete Bordeaux input:", len(fields))
    for row in square_rows:
        print(row["discriminant"], row["polynomial"], "h =", row["ordinary_class_number"])
    print("E2 square-ramified quartic certificate: all checks passed")


if __name__ == "__main__":
    main()
