"""Apply exact all-component lattice certificates to an enumerated field list.

Run the candidate enumerator first, then:
  sage -python e2_indecomposable_search.sage --input e2_candidates.json --output e2_lattices.json

A missing prime below --prime-bound produces an explicit unresolved record
and a nonzero exit status.  It is never silently treated as a proof.
The four named arithmetic exceptions are checked up to field isomorphism
and referred to the separate dyadic/space computations, not excluded here.
"""

from sage.all import *
from sage.version import version as SAGE_VERSION
from pathlib import Path
import argparse
import json
import sys

proof.all(True)

# .sage files here deliberately use ordinary Python syntax; exec avoids
# Sage preparser and path-dependent load side effects in the module import.
module_scope={"__name__":"case11_lattice_module"}
exec(compile((Path(__file__).resolve().parent/"e2_allclass_balanced.sage").read_text(),
             "e2_allclass_balanced.sage","exec"),module_scope)
find_certificate=module_scope["find_certificate"]

EXCEPTION_POLYNOMIALS={
    (3,49):[1,-2,-1,1],
    (4,1125):[1,4,-4,-1,1],
    (4,5125):[11,7,-6,-2,1],
    (4,6125):[11,9,-9,-1,1],
}


def is_arithmetic_exception(F):
    key=(int(F.degree()),int(F.discriminant()))
    if key not in EXCEPTION_POLYNOMIALS:
        return False
    Qx=PolynomialRing(QQ,"x")
    other=NumberField(Qx(EXCEPTION_POLYNOMIALS[key]),"b")
    assert other.discriminant()==F.discriminant()
    return bool(F.is_isomorphic(other))


def run_manifest(payload, prime_bound):
    if not payload.get("enumeration_complete",False):
        raise ValueError("incomplete enumeration: resume/re-run it before lattice certification")
    Qx=PolynomialRing(QQ,"x")
    results=[]
    for source in payload["fields"]:
        if source["status"]!="requires_lattice_or_exception_certificate":
            continue
        F=NumberField(Qx(source["polynomial"]),"a")
        assert F.discriminant()==source["discriminant"]
        assert F.class_group(proof=True).order()==source["ordinary_class_number"]
        row=dict(source)
        if is_arithmetic_exception(F):
            row["status"]="external_arithmetic_or_space_certificate_required"
        else:
            certificate,attempts=find_certificate(F,prime_bound)
            row["certificate"]=certificate
            row["attempts"]=attempts
            row["status"]="excluded_by_lattice_certificate" if certificate else "unresolved_search"
        results.append(row)
        print(F.degree(),F.discriminant(),row["status"],flush=True)
        yield results


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",default="e2_candidates.json")
    parser.add_argument("--output",default="e2_lattices.json")
    parser.add_argument("--prime-bound",type=int,default=43)
    args=parser.parse_args()
    source=json.loads(Path(args.input).read_text())
    output={"schema":"case11-lattice-search-v1","sage_version":SAGE_VERSION,
            "input_manifest":args.input,"prime_bound":args.prime_bound,
            "all_required_degrees":source.get("covers_all_required_degrees",False),
            "search_complete":False,"fields":[]}
    for rows in run_manifest(source,args.prime_bound):
        output["fields"]=rows
        Path(args.output).write_text(json.dumps(output,indent=2))
    output["search_complete"]=True
    unresolved=[r for r in output["fields"] if r["status"]=="unresolved_search"]
    output["unresolved_fields"]=[{"degree":r["degree"],"discriminant":r["discriminant"],
        "polynomial":r["polynomial"]} for r in unresolved]
    output["all_nonexception_candidates_certified"]=not unresolved and output["all_required_degrees"]
    Path(args.output).write_text(json.dumps(output,indent=2))
    print("Nonexception lattice certification complete:",output["all_nonexception_candidates_certified"])
    if unresolved or not output["all_required_degrees"]:
        sys.exit(2)
