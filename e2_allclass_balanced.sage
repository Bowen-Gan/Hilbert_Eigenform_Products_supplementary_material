"""Exact positive decomposition certificates, in every narrow class.

Run: sage -python e2_allclass_balanced.sage --output d3969_certificate.json

For any integral ideal A, choose the positive integer m generating A cap Z
and set t=(m) A^-1.  Then A=(m)t^-1 and t is integral.  Testing whether
m decomposes in t does not depend on this choice: simultaneous positive
scaling of m and t preserves every decomposition.  No unit search, narrow
principality assumption, or truncated short-vector search is required.

The only use of LLL is an invertible integral change of basis.  Its
unimodularity is checked exactly; its quality is irrelevant to correctness.
Every integer point in an exact trace-dual box is tested with algebraic
real embeddings.  Exhaustion of that box certifies indecomposability.
"""

from sage.all import *
from sage.version import version as SAGE_VERSION
from itertools import product
import argparse
import json

proof.all(True)


def field_element_data(a):
    return [str(c) for c in a.list()]


def matrix_data(M):
    return [[str(c) for c in row] for row in M.rows()]


def trace_basis(ideal, reduce_basis=True):
    F = ideal.number_field()
    basis = list(ideal.basis())
    n = F.degree()
    gram = matrix(ZZ, [[(u*v).trace() for v in basis] for u in basis])
    assert gram.is_positive_definite()
    U = identity_matrix(ZZ, n)
    if reduce_basis:
        U = matrix(ZZ, pari(gram).qflllgram())
        assert abs(U.det()) == 1
        basis = [sum((basis[i]*U[i,j] for i in range(n)), F(0))
                 for j in range(n)]
        gram = U.transpose()*gram*U
        assert gram == matrix(ZZ, [[(u*v).trace() for v in basis] for u in basis])
    assert F.ideal(basis) == ideal
    return basis, gram, U


def coordinate_box(mu, basis, gram):
    F = mu.parent()
    # Outward intervals suffice for a containing box.  Computing sums of
    # unrelated exact conjugates in AA can unnecessarily construct a large
    # normal closure; rational interval endpoints avoid that expense.
    embeddings = F.embeddings(RealIntervalField(160))
    assert len(embeddings) == F.degree()
    mu_values = [s(mu) for s in embeddings]
    assert all(v > 0 for v in mu_values)
    inv = gram.inverse()
    dual = [sum((inv[j,k]*basis[k] for k in range(F.degree())), F(0))
            for j in range(F.degree())]
    assert matrix(QQ, [[(u*v).trace() for v in dual] for u in basis]).is_one()
    bounds = []
    for bstar in dual:
        values = [m*s(bstar) for m,s in zip(mu_values,embeddings)]
        if all(v.lower()>0 for v in values):
            lo,hi=QQ(0),QQ((mu*bstar).trace())
        elif all(v.upper()<0 for v in values):
            lo,hi=QQ((mu*bstar).trace()),QQ(0)
        else:
            lo = sum((min(QQ(0),QQ(v.lower())) for v in values),QQ(0))
            hi = sum((max(QQ(0),QQ(v.upper())) for v in values),QQ(0))
        # An enlarged closed box is used so that the transcript is simple.
        # 0 << x << mu implies lo < Tr(x bstar) < hi, hence inclusion.
        bounds.append((ZZ(lo.floor()), ZZ(hi.ceil())))
    return bounds, embeddings


def certify_exponent(ideal, include_points=False):
    F = ideal.number_field()
    mu = F(ideal.smallest_integer())
    assert mu > 0 and mu in ideal
    component = F.ideal(mu)*ideal**(-1)
    assert component.is_integral() and mu in component
    assert F.ideal(mu)*component**(-1) == ideal
    basis, gram, U = trace_basis(component)
    bounds, embeddings = coordinate_box(mu, basis, gram)
    expected = prod(b-a+1 for a,b in bounds)
    checked = 0
    witness = None
    points = []
    for coords in product(*(range(int(a), int(b)+1) for a,b in bounds)):
        checked += 1
        z = sum((ZZ(c)*b for c,b in zip(coords, basis)), F(0))
        def exact_signs(element):
            if element==0:
                return [0]*F.degree()
            result=[]
            for s in embeddings:
                value=s(element)
                if value>0:
                    result.append(1)
                elif value<0:
                    result.append(-1)
                else:
                    # An ambiguous interval cannot certify a sign.  Exact
                    # algebraic signs are used for the rare close-to-zero case.
                    return [int(sign(e(element))) for e in F.embeddings(AA)]
            return result
        signs_z = exact_signs(z)
        signs_rest = exact_signs(mu-z)
        if include_points:
            points.append({"coordinates": list(coords), "signs_x": signs_z,
                           "signs_mu_minus_x": signs_rest})
        if all(v == 1 for v in signs_z) and all(v == 1 for v in signs_rest):
            witness = field_element_data(z)
            break
    status = "indecomposable" if witness is None else "decomposable"
    if witness is None:
        assert checked == expected
    return {"status": status, "mu": field_element_data(mu),
            "ideal_basis": [field_element_data(b) for b in ideal.basis()],
            "component_basis": [field_element_data(b) for b in basis],
            "basis_change": matrix_data(U), "trace_gram": matrix_data(gram),
            "coordinate_bounds": [[int(a),int(b)] for a,b in bounds],
            "box_size": int(expected), "points_checked": checked,
            "decomposition": witness, "points": points}


def certify_prime(P, include_points=False):
    assert P.is_prime()
    data = []
    for exponent in (1,2):
        row = certify_exponent(P**exponent, include_points=include_points)
        row["exponent"] = exponent
        data.append(row)
        if row["status"] != "indecomposable":
            return None, data
    return {"prime_basis": [field_element_data(b) for b in P.basis()],
            "norm": int(P.norm()), "exponents": data}, data


def find_certificate(F, rational_bound=43):
    prime_rows = []
    for p in prime_range(3, rational_bound+1):
        for P,e in F.ideal(p).factor():
            prime_rows.append((ZZ(P.norm()), ZZ(p), str(P), P, ZZ(e)))
    prime_rows.sort(key=lambda r: r[:3])
    attempts = []
    for normP,p,_,P,e in prime_rows:
        certificate, data = certify_prime(P)
        attempts.append({"p": int(p), "norm": int(normP), "e": int(e),
                         "exponents": data})
        if certificate is not None:
            certificate.update({"rational_prime": int(p), "ramification_index": int(e)})
            return certificate, attempts
    # This is a bounded search failure, never an exclusion or proof of absence.
    return None, attempts


def d3969_certificate():
    x = polygen(QQ)
    F = NumberField(x**3-21*x-35, "a")
    a = F.gen()
    assert F.discriminant() == 3969
    assert tuple(F.class_group(proof=True).invariants()) == (3,)
    assert tuple(F.narrow_class_group().invariants()) == (3,)
    P = F.ideal(3,a+1)
    assert P.is_prime() and P.norm() == 3 and P**3 == F.ideal(3)
    certificate, data = certify_prime(P, include_points=True)
    assert certificate is not None
    assert all(r["mu"] == ["3", "0", "0"] or r["mu"] == ["3"] for r in data)
    return {"polynomial": [-35,-21,0,1], "discriminant":3969,
            "class_group":[3], "narrow_class_group":[3], "certificate":certificate}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="d3969_certificate.json")
    args = parser.parse_args()
    result = d3969_certificate()
    result["sage_version"] = SAGE_VERSION
    with open(args.output,"w") as handle:
        json.dump(result, handle, indent=2)
    print("D=3969: both exponents indecomposable; output:", args.output)
