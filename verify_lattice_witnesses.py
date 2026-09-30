#!/usr/bin/env python3
r"""Verify the fifteen stored Case 11 lattice certificates independently.

Standard library only: rational polynomial arithmetic, Sturm root isolation,
integer lattice indices, and outward rational intervals. Every integer point
in each containing trace-dual box is checked; no floating-point sign decision
and no Sage installation is used. A failed check exits nonzero.

Examples (run from the supplement root):
    python verify_lattice_witnesses.py
    python verify_lattice_witnesses.py \
        --data-dir rerun --output lattice_verification.json

External inputs remain explicit. Except for discriminant 3969, identification
of the supplied integer orders with the maximal orders relies on the reported
field discriminants. The verifier checks the ring and lattice arithmetic and
that each order's discriminant equals its input field discriminant; it does
not independently enumerate number fields or recompute all discriminants or
class groups. At discriminant 3969, Dedekind's index criterion independently
establishes maximality. The Hilbert-modular coefficient comparison and the
analytic dimension criterion are theoretical inputs, not checked here.

Only the fifteen lattice exclusions are certified. The six other candidates
are reported as requiring their separate arithmetic or space arguments; they
are never promoted to lattice exclusions.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from math import prod, isqrt
import argparse
import hashlib
import json
import sys

# Assertions are part of the verifier, so optimized Python must be rejected.
if not __debug__:
    print("ERROR: optimized Python disables verification assertions; rerun without -O/-OO or PYTHONOPTIMIZE.", file=sys.stderr)
    raise SystemExit(2)

def trim(a):
    while len(a)>1 and a[-1]==0: a.pop()
    return a

def divrem(a,b):
    a=list(a); q=[Q(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and any(a):
        k=len(a)-len(b); c=a[-1]/b[-1];q[k]=c
        for j,x in enumerate(b):a[j+k]-=c*x
        trim(a)
    return trim(q),trim(a)

def val(a,x):
    y=Q(0)
    for c in reversed(a):y=y*x+c
    return y

def sign(x):return (x>0)-(x<0)

def intervals(f):
    st=[f,[j*f[j] for j in range(1,len(f))]]
    while len(st[-1])>1:
        r=divrem(st[-2],st[-1])[1]; assert any(r)
        st.append([-x for x in r])
    def variation(x):
        s=[sign(val(p,x)) for p in st];s=[x for x in s if x]
        return sum(a!=b for a,b in zip(s,s[1:]))
    M=1+max(abs(x) for x in f[:-1]);todo=[(-M,M)];ans=[]
    while todo:
        a,b=todo.pop();count=variation(a)-variation(b)
        if not count:continue
        if count==1 and b-a<Q(1,2**45):ans.append((a,b));continue
        m=(a+b)/2;assert val(f,m)!=0
        todo.extend([(a,m),(m,b)])
    assert len(ans)==len(f)-1
    return sorted(ans)

def iadd(a,b):return (a[0]+b[0],a[1]+b[1])
def imul(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)
def iscale(c,a):return (c*a[0],c*a[1]) if c>=0 else (c*a[1],c*a[0])
def ival(p,r):
    y=(Q(0),Q(0))
    for c in reversed(p):y=iadd(imul(y,r),(c,c))
    return y
def isign(a):
    if a[0]>0:return 1
    if a[1]<0:return -1
    raise AssertionError(('undecided interval',a))

def inverse(a):
    n=len(a);m=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(a)]
    for k in range(n):
        i=next(i for i in range(k,n) if m[i][k]);m[k],m[i]=m[i],m[k]
        v=m[k][k];m[k]=[x/v for x in m[k]]
        for i in range(n):
            if i!=k:
                v=m[i][k];m[i]=[x-v*y for x,y in zip(m[i],m[k])]
    return [r[n:] for r in m]

def det(a):
    a=[list(r) for r in a];v=Q(1);n=len(a)
    for k in range(n):
        if not any(a[i][k] for i in range(k,n)):return Q(0)
        i=next(i for i in range(k,n) if a[i][k])
        if i!=k:a[k],a[i]=a[i],a[k];v=-v
        t=a[k][k];v*=t
        for i in range(k+1,n):
            c=a[i][k]/t
            for j in range(k+1,n):a[i][j]-=c*a[k][j]
    return v

def rowmul(v,a):return [sum(x*a[i][j] for i,x in enumerate(v)) for j in range(len(a))]
def lattice_index(rows,n):
    a=[]
    for row in rows:
        assert all(x.denominator==1 for x in row)
        a.append([int(x) for x in row])
    for k in range(n):
        assert any(a[i][k] for i in range(k,len(a)))
        while True:
            i=min((i for i in range(k,len(a)) if a[i][k]),key=lambda i:abs(a[i][k]))
            a[k],a[i]=a[i],a[k]
            for i in range(k+1,len(a)):
                q=a[i][k]//a[k][k]
                a[i]=[x-q*y for x,y in zip(a[i],a[k])]
            if all(a[i][k]==0 for i in range(k+1,len(a))):break
    assert all(not any(r) for r in a[n:])
    return abs(prod(a[k][k] for k in range(n)))

def matrix(data,n):return [[Q(x) for x in r]+[Q(0)]*(n-len(r)) for r in data]

class Field:
    def __init__(self,data):
        self.n=data['degree'];self.f=list(map(Q,data['polynomial']));n=self.n
        assert len(self.f)==n+1 and self.f[-1]==1
        # Every selected defining polynomial is irreducible modulo 2.
        def remainder2(f,g):
            f=list(f)
            while len(f)>=len(g):
                if f[-1]:
                    k=len(f)-len(g)
                    for j,c in enumerate(g):f[k+j]^=c
                trim(f)
            return f
        ff=[int(x)%2 for x in self.f]
        for d in range(1,n//2+1):
            for bits in product([0,1],repeat=d):
                assert any(remainder2(ff,list(bits)+[1]))
        self.roots=intervals(self.f)
        self.traces=[Q(n)]
        for k in range(1,n):
            self.traces.append(-sum(self.f[n-j]*self.traces[k-j] for j in range(1,k))-k*self.f[n-k])
        dyad=data['dyadic_factorization'];assert len(dyad)==1
        assert dyad[0]['norm']==2**n and dyad[0]['ramification_index']==1
        self.O=[[x/2 for x in r] for r in matrix(dyad[0]['basis'],n)]
        self.Oinv=inverse(self.O)
        assert det(self.gram(self.O))==data['discriminant']
        if data['discriminant']==3969:
            # Only 3 and 7 can divide the power-basis index. Modulo 3 the
            # polynomial is (x+1)^3; modulo 7 it is x^3. In each case
            # (f-g^3)/p is nonzero at the unique repeated root, so
            # Dedekind's index criterion excludes p from the index.
            assert self.f==list(map(Q,[-35,-21,0,1]))
            assert 3969==3**4*7**2
            for p,gcube,root in [(3,[1,3,3,1],-1),(7,[0,0,0,1],0)]:
                quotient=[(a-b)/p for a,b in zip(self.f,gcube)]
                assert all(x.denominator==1 for x in quotient)
                assert val(quotient,Q(root))%p != 0
        assert all(x.denominator==1 for x in rowmul([Q(1)]+[Q(0)]*(n-1),self.Oinv))
        for a in self.O:
            for b in self.O:
                assert all(x.denominator==1 for x in rowmul(self.mul(a,b),self.Oinv))
    def mul(self,a,b):
        n=self.n;c=[Q(0)]*(2*n-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):c[i+j]+=x*y
        r=divrem(c,self.f)[1];return r+[Q(0)]*(n-len(r))
    def trace(self,a):return sum(x*y for x,y in zip(a,self.traces))
    def gram(self,b):return [[self.trace(self.mul(x,y)) for y in b] for x in b]
    def check_ideal(self,I):
        inv=inverse(I)
        for b in I:assert all(x.denominator==1 for x in rowmul(b,self.Oinv))
        for a in self.O:
            for b in I:assert all(x.denominator==1 for x in rowmul(self.mul(a,b),inv))
        return inv
    def product_equals(self,A,B,C):
        inv=inverse(C)
        assert lattice_index([rowmul(self.mul(a,b),inv) for a in A for b in B],self.n)==1

def verify(data,sign_data=None):
    F=Field(data);n=F.n;cert=data['certificate'];P=matrix(cert['prime_basis'],n)
    assert [r['exponent'] for r in cert['exponents']]==[1,2]
    F.check_ideal(P);p=cert['norm']
    assert abs(det(P)/det(F.O))==p
    assert all(p%d for d in range(2,isqrt(p)+1))
    total=0;box_rows=[]
    for row in cert['exponents']:
        exponent=row['exponent'];I=matrix(row['ideal_basis'],n);B=matrix(row['component_basis'],n)
        F.check_ideal(I);F.check_ideal(B)
        if exponent==1:assert lattice_index([rowmul(x,inverse(P)) for x in I],n)==1
        else:
            F.product_equals(P,P,I)
            if data['discriminant']==3969:
                F.product_equals(P,I,[[3*x for x in b] for b in F.O])
        mu=list(map(Q,row['mu']));mu += [Q(0)]*(n-len(mu));assert mu[0]>0 and not any(mu[1:])
        assert all(x.denominator==1 for x in rowmul(mu,inverse(I)))
        assert all(x.denominator==1 for x in rowmul(mu,inverse(B)))
        F.product_equals(B,I,[[mu[0]*x for x in b] for b in F.O])
        U=matrix(row['basis_change'],n);assert abs(det(U))==1 and all(x.denominator==1 for r in U for x in r)
        gram=F.gram(B);assert gram==matrix(row['trace_gram'],n)
        inv=inverse(gram);dual=[rowmul(r,B) for r in inv]
        for i in range(n):
            for j in range(n):assert F.trace(F.mul(B[i],dual[j]))==(i==j)
        exact_bounds=[]
        for d in dual:
            vals=[iscale(mu[0],ival(d,r)) for r in F.roots]
            signs=[isign(v) for v in vals]
            if all(s==1 for s in signs):lo=Q(0);hi=mu[0]*F.trace(d)
            elif all(s==-1 for s in signs):lo=mu[0]*F.trace(d);hi=Q(0)
            else:
                lo=sum(min(Q(0),v[0]) for v in vals)
                hi=sum(max(Q(0),v[1]) for v in vals)
            exact_bounds.append([lo.numerator//lo.denominator,-((-hi.numerator)//hi.denominator)])
        assert exact_bounds==row['coordinate_bounds'],(data['discriminant'],exponent,exact_bounds,row['coordinate_bounds'])
        bounds=row['coordinate_bounds'];count=prod(b-a+1 for a,b in bounds)
        assert count==row['box_size']==row['points_checked']
        embeds=[[ival(b,r) for r in F.roots] for b in B]
        points=[]
        for coords in product(*(range(a,b+1) for a,b in bounds)):
            element=rowmul(coords,B);rest=[x-y for x,y in zip(mu,element)]
            vals=[(Q(0),Q(0)) for _ in range(n)]
            for c,e in zip(coords,embeds):
                vals=[iadd(v,iscale(c,w)) for v,w in zip(vals,e)]
            sx=[0]*n if not any(element) else [isign(v) for v in vals]
            sy=[0]*n if not any(rest) else [isign((mu[0]-v[1],mu[0]-v[0])) for v in vals]
            assert not (all(s==1 for s in sx) and all(s==1 for s in sy))
            points.append({'coordinates':list(coords),'signs_x':sx,'signs_mu_minus_x':sy})
        if row['points']:assert points==row['points']
        if sign_data is not None:assert points==sign_data['certificate']['exponents'][exponent-1]['points']
        total+=count;box_rows.append(count)
    print('PASS',data['degree'],data['discriminant'],'boxes',box_rows,flush=True)
    return {
        'degree': data['degree'],
        'discriminant': data['discriminant'],
        'polynomial': data['polynomial'],
        'prime_norm': int(p),
        'exponent_boxes': len(box_rows),
        'box_sizes': box_rows,
        'checked_points': total,
        'indecomposable_exponents': [1, 2],
        'maximality_basis': ('independent Dedekind index criterion' if data['discriminant']==3969
                            else 'order discriminant matches externally supplied field discriminant'),
        'passed': True,
    }

EXPECTED_LATTICE_FIELDS = {
    (3,81),(3,257),(3,321),(3,697),(3,1257),(3,1489),(3,3969),
    (4,2525),(4,4205),(4,4525),(4,8069),(4,16317),
    (5,14641),(5,38569),(6,371293),
}
EXPECTED_SEPARATE_FIELDS = {(3,49),(3,169),(3,361),(4,725),(4,1125),(4,5125)}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(data_dir):
    lattice_path=data_dir/'e2_lattices.json'
    special_path=data_dir/'d3969_certificate.json'
    payload=json.loads(lattice_path.read_text())
    special=json.loads(special_path.read_text())
    assert payload['search_complete'] is True
    all_rows=payload['fields']
    assert len(all_rows)==21
    rows=[r for r in all_rows if r.get('certificate')]
    separate=[r for r in all_rows if not r.get('certificate')]
    assert len(rows)==15
    assert {(r['degree'],r['discriminant']) for r in rows}==EXPECTED_LATTICE_FIELDS
    assert {(r['degree'],r['discriminant']) for r in separate}==EXPECTED_SEPARATE_FIELDS
    assert all(r['status']=='excluded_by_lattice_certificate' for r in rows)
    unresolved=[r for r in all_rows if r['status']=='unresolved_search']
    assert {(r['degree'],r['discriminant']) for r in unresolved}=={(3,169),(3,361),(4,725)}
    # Preserve the search's honest unresolved flag; separate arguments must
    # establish these fields' exclusions and are outside this verifier.
    expected_flag=not unresolved and bool(payload['all_required_degrees'])
    assert payload['all_nonexception_candidates_certified']==expected_flag
    assert payload['all_nonexception_candidates_certified'] is False
    assert special['discriminant']==3969 and special['polynomial']==[-35,-21,0,1]
    field_results=[verify(r,special if r['discriminant']==3969 else None) for r in rows]
    total=sum(r['checked_points'] for r in field_results)
    boxes=sum(r['exponent_boxes'] for r in field_results)
    assert boxes==30 and total==8776
    sign_records=sum(len(r['points']) for r in special['certificate']['exponents'])
    assert sign_records==36
    return {
        'schema': 'case11-independent-lattice-verification-v1',
        'status': 'passed', 'completed': True, 'passed': True,
        'scope': 'fifteen supplied lattice exclusions; thirty containing boxes',
        'python_version': sys.version.split()[0],
        'arithmetic': 'integer/rational arithmetic; Sturm-isolated roots; outward rational intervals',
        'source_hashes': {p.name: sha256(p) for p in [lattice_path,special_path]},
        'verifier_sha256': sha256(Path(__file__).resolve()),
        'field_count': len(field_results), 'exponent_box_count': boxes,
        'checked_point_count': total,
        'D3969_sign_record_count': sign_records,
        'D3969_maximal_order_independently_checked': True,
        'source_lattice_search_all_nonexception_candidates_certified':
            payload['all_nonexception_candidates_certified'],
        'fields': field_results,
        'separate_arguments_required': [
            {'degree': r['degree'], 'discriminant': r['discriminant'],
             'source_status': r['status'], 'lattice_exclusion_verified_here': False}
            for r in separate],
        'external_dependencies': [
            'Complete enumeration of the original 771 fields is not reproduced here.',
            'Reported field discriminants establish maximal-order identification outside D=3969; class groups are not recomputed.',
            'The small-norm/indecomposable Fourier coefficient comparison and Hecke recurrences are mathematical inputs.',
            'The analytic all-component pairing and dimension criterion are not verified by this computation.',
            'The six remaining candidates require separate arithmetic, dyadic, or cusp-space arguments.',
        ],
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parent/'rerun',
                        help='directory containing e2_lattices.json and d3969_certificate.json (default: repository rerun directory)')
    parser.add_argument('--output', type=Path, help='optional JSON result path')
    args=parser.parse_args()
    try:
        result=run(args.data_dir.resolve())
    except Exception as exc:
        result={'schema':'case11-independent-lattice-verification-v1',
                'status':'failed','completed':False,'passed':False,
                'error_type':type(exc).__name__,'error':str(exc)}
        if args.output:
            args.output.write_text(json.dumps(result,indent=2)+'\n')
        print('FAIL: '+type(exc).__name__+': '+str(exc),file=sys.stderr)
        return 1
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('ALL_PASS: fields=%d exponent_boxes=%d checked_points=%d D3969_sign_records=%d D3969_maximal_order=True' % (
        result['field_count'], result['exponent_box_count'], result['checked_point_count'],result['D3969_sign_record_count']))
    print('SCOPE: six separate candidates are not certified as lattice exclusions; external dependencies remain explicit.')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
