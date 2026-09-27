"""Complete finite input for the Case 11 lattice calculations.

Run: sage -python e2_candidate_enumeration.sage --output e2_candidates.json

This enumerates fields, not discriminants, and retains each defining
polynomial.  No list of 'residual discriminants' is assumed.  The routine
enumerate_totallyreal_fields_all includes imprimitive fields.

Scope: apply the manuscript's general even-class-number exclusion first.
If h is odd and n*h <= 14, then h>1 is possible only for n=3,4.
Thus n=3,4 require the general alpha bound; n=5,6 only require the
sharper h=1 bound.  Degrees >=7 are handled in the high-degree proof.

The integer endpoints below are certified using outward interval arithmetic.
They are intentionally rounded upwards; every enumerated field is recorded,
including those subsequently excluded.  A status describes the manuscript
criterion used, not an independent proof of that criterion.
"""

from sage.all import *
from sage.version import version as SAGE_VERSION
from sage.rings.number_field.totallyreal_rel import enumerate_totallyreal_fields_all
import argparse
import json

proof.all(True)

ENUMERATION_BOUNDS = {3:4218, 4:68133, 5:209508, 6:2429064}
ENUMERATION_SOURCE = "https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/totallyreal_rel.html"


def exact_endpoint_checks():
    R = RealIntervalField(256)
    pi = R.pi()
    for n,B in ENUMERATION_BOUNDS.items():
        base = 2*pi**4/3 if n in (3,4) else 4*pi**2
        assert R(B)**3 > base**(2*n)
    return {str(n): B for n,B in ENUMERATION_BOUNDS.items()}


def local_record(F,p):
    return [{"ramification_index": int(e), "residue_degree": int(P.residue_class_degree()),
             "norm": int(P.norm()), "basis": [[str(c) for c in b.list()] for b in P.basis()]}
            for P,e in F.ideal(p).factor()]


def reciprocal_zeta_interval(F, cutoff=199):
    """Rigorous |2^n/zeta_F(-1)| interval from finite Euler factors.

    P_B < zeta_F(2) < P_B*((B+1)/B)^n.  The upper tail follows by
    bounding every rational-prime Euler factor by (1-p^-2)^-n, then
    using the telescoping product over ALL integers m>B.  No floating
    special-value reconstruction is used.
    """
    n=F.degree()
    P=QQ(1)
    factors=[]
    for p in prime_range(2,cutoff+1):
        degrees=[ZZ(prime.residue_class_degree()) for prime,e in F.ideal(p).factor()]
        for f in degrees:
            P /= 1-QQ(p)**(-2*f)
        factors.append([int(p),[int(f) for f in degrees]])
    tail=(QQ(cutoff+1)/cutoff)**n
    R=RealIntervalField(256)
    D=R(F.discriminant())
    C=(4*R.pi()**2)**n/(D*D.sqrt())
    lower=(C/(R(P)*R(tail))).lower()
    upper=(C/R(P)).upper()
    lo=ZZ(lower.ceil())
    hi=ZZ(upper.floor())
    return {"cutoff":cutoff,"factor_degrees":factors,
            "euler_product":str(P),"tail_bound":str(tail),
            "alpha_absolute_lower":str(lower),"alpha_absolute_upper":str(upper),
            "possible_integer_min":int(lo),"possible_integer_max":int(hi),
            "excludes_integrality":lo>hi}


def classify_field(F):
    n = F.degree()
    D = ZZ(F.discriminant())
    h = ZZ(F.class_group(proof=True).order())
    row = {"degree":int(n), "discriminant":int(D),
           "polynomial":[int(c) for c in F.defining_polynomial().list()],
           "ordinary_class_number":int(h)}
    if h % 2 == 0:
        row["status"] = "even_class_number_argument"
        return row
    if n*h > 14:
        row["status"] = "general_degree_class_number_bound"
        return row
    R = RealIntervalField(256)
    if h == 1 and R(D)**3 > (4*R.pi()**2)**(2*n):
        row["status"] = "trivial_character_discriminant_bound"
        return row
    row["dyadic_factorization"] = local_record(F,2)
    facts = row["dyadic_factorization"]
    if not (len(facts)==1 and facts[0]["ramification_index"]==1
            and facts[0]["residue_degree"]==n):
        row["status"] = "noninert_two_argument"
        return row
    for p in prime_range(3, ZZ(2**n).isqrt()+1):
        for P,e in F.ideal(p).factor():
            if P.norm()**2 <= 2**n and P**2 != F.ideal(2):
                row["status"] = "small_prime_argument"
                row["exclusion_prime"] = {"p":int(p), "norm":int(P.norm()),
                    "basis":[[str(c) for c in b.list()] for b in P.basis()]}
                return row
    if h==1:
        row["alpha_interval"]=reciprocal_zeta_interval(F)
        if row["alpha_interval"]["excludes_integrality"]:
            row["status"]="reciprocal_constant_not_integral"
            return row
    row["narrow_class_number"] = int(F.narrow_class_group().order())
    row["status"] = "requires_lattice_or_exception_certificate"
    return row


def enumerate_candidates(degrees=(3,4,5,6), checkpoint=None):
    exact_endpoint_checks()
    results = []
    for n in degrees:
        if n not in ENUMERATION_BOUNDS:
            raise ValueError("only degrees 3,4,5,6 are in this finite reduction")
        print("Enumerating degree",n,"through",ENUMERATION_BOUNDS[n],flush=True)
        fields = enumerate_totallyreal_fields_all(n, ENUMERATION_BOUNDS[n],
                                                  return_pari_objects=False)
        seen = set()
        for D,f in fields:
            F = NumberField(f,"a")
            assert F.degree()==n and F.signature()==(n,0)
            assert F.discriminant()==D
            key = tuple(ZZ(c) for c in F.defining_polynomial().list())
            assert key not in seen
            seen.add(key)
            row = classify_field(F)
            results.append(row)
            print(n,D,row["ordinary_class_number"],row["status"],flush=True)
            if checkpoint is not None:
                checkpoint(results, list(degrees), n, False)
    return results


def write_manifest(path, rows, degrees, current_degree=None, complete=False):
    payload = {"schema":"case11-field-enumeration-v1", "sage_version":SAGE_VERSION,
               "enumeration_function":"enumerate_totallyreal_fields_all",
               "source":ENUMERATION_SOURCE, "endpoint_checks":exact_endpoint_checks(),
               "requested_degrees":degrees, "current_degree":current_degree,
               "enumeration_complete":complete,
               "covers_all_required_degrees":complete and set(degrees)=={3,4,5,6},
               "fields":rows}
    with open(path,"w") as handle:
        json.dump(payload,handle,indent=2)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",default="e2_candidates.json")
    parser.add_argument("--degrees",nargs="+",type=int,default=[3,4,5,6])
    args=parser.parse_args()
    rows=enumerate_candidates(args.degrees,
        checkpoint=lambda rows,degrees,n,complete:write_manifest(args.output,rows,degrees,n,complete))
    write_manifest(args.output,rows,args.degrees,complete=True)
    remaining=[r for r in rows if r["status"]=="requires_lattice_or_exception_certificate"]
    print("Enumeration complete; lattice/exception inputs:",len(remaining))
    print("Odd class number >1:",[(r["degree"],r["discriminant"],r["ordinary_class_number"])
                                      for r in remaining if r["ordinary_class_number"]>1])
