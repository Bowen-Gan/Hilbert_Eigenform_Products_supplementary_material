#!/usr/bin/env python3
"""Independent rational verification of the six dyadic weight cutoffs.

Python standard library only. Reconstructs trace matrices, exhaustive boxes,
positive decompositions, principal-ideal factors, and the Euler lower products
from the defining polynomials. Decisions use integers and rational intervals.

External inputs: the supplied polynomials have the stated FIELD discriminants;
the six exact/upper reciprocal constants and the theoretical coefficient bound.
Equality of the power-basis discriminant with the field discriminant then proves
that the power basis is maximal, so Dedekind factorization applies at every p.
For 49,169,361,1125 the stated exact alpha values remain special-value inputs.
For 725,5125 the alpha upper bounds are independently derived below.
This program does not verify Rankin--Selberg or Ramanujan theorems.
"""
import argparse
import hashlib
import json
import sys
from fractions import Fraction as Q
from itertools import product
from math import isqrt, prod
from pathlib import Path

FIELDS = {
    49: ([1, -2, -1, 1], Q(168), 7, 34),
    169: ([-1, -4, -1, 1], Q(24), 4, 39),
    361: ([7, -6, -1, 1], Q(8), 3, 34),
    725: ([1, 1, -3, -1, 1], Q(121), 4, 138),
    1125: ([1, 4, -4, -1, 1], Q(60), 4, 82),
    5125: ([11, 7, -6, -2, 1], Q(61, 10), 2, 149),
}


def trim(f):
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f


def divrem(f, g, p=None):
    f = [Q(x) for x in f] if p is None else [int(x) % p for x in f]
    g = [Q(x) for x in g] if p is None else [int(x) % p for x in g]
    q = [0] * max(1, len(f) - len(g) + 1)
    while len(f) >= len(g) and any(f):
        j = len(f) - len(g)
        c = f[-1] / g[-1] if p is None else f[-1] * pow(g[-1], -1, p) % p
        q[j] = c
        for i, a in enumerate(g):
            f[i+j] -= c*a
            if p is not None:
                f[i+j] %= p
        trim(f)
    return trim(q), trim(f)


def mul(a, b, f, p=None):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return divrem(c, f, p)[1]


def gcdpoly(a, b, p):
    a, b = trim(a[:]), trim(b[:])
    while any(b):
        a, b = b, divrem(a, b, p)[1]
    c = pow(a[-1], -1, p)
    return [(int(x)*c) % p for x in a]


def powpoly(a, e, f, p):
    v = [1]
    while e:
        if e % 2:
            v = mul(v, a, f, p)
        a = mul(a, a, f, p)
        e //= 2
    return v


def factor_norms(f, p):
    """Complete residue degrees and ramification indices for degree <=4."""
    r = [x % p for x in f]
    out = []
    for a in range(p):
        g = [(-a) % p, 1]
        e = 0
        while len(r) > 1 and not any(divrem(r, g, p)[1]):
            r = divrem(r, g, p)[0]
            e += 1
        if e:
            out.append((p, e))
    d = len(r)-1
    if d in (2, 3):
        # A quadratic or cubic without a linear factor is irreducible.
        out.append((p**d, 1))
    elif d == 4:
        derivative = [(i*r[i]) % p for i in range(1, len(r))]
        repeated = gcdpoly(r, derivative, p) if any(derivative) else r
        if len(repeated) > 1:
            # A root-free reducible repeated quartic is a quadratic square.
            assert len(repeated) in (3, 5)
            out.append((p**2, 2))
        else:
            xpp = powpoly([0, 1], p*p, r, p)
            xpp += [0] * max(0, 2-len(xpp))
            xpp[1] = (xpp[1]-1) % p
            d2 = len(gcdpoly(r, xpp, p))-1
            assert d2 in (0, 4)
            out.extend([(p**2, 1), (p**2, 1)] if d2 == 4 else [(p**4, 1)])
    else:
        assert d == 0
    assert prod(q**e for q, e in out) == p**(len(f)-1)
    return sorted(out)


def det(a):
    a = [[Q(x) for x in r] for r in a]
    z = Q(1)
    for k in range(len(a)):
        j = next((j for j in range(k, len(a)) if a[j][k]), None)
        if j is None:
            return Q(0)
        if j != k:
            a[j], a[k] = a[k], a[j]
            z = -z
        t = a[k][k]
        z *= t
        for j in range(k+1, len(a)):
            v = a[j][k]/t
            for i in range(k+1, len(a)):
                a[j][i] -= v*a[k][i]
    return z


def inverse(a):
    n = len(a)
    a = [[Q(x) for x in a[i]] + [Q(i == j) for j in range(n)] for i in range(n)]
    for k in range(n):
        j = next(j for j in range(k, n) if a[j][k])
        a[j], a[k] = a[k], a[j]
        t = a[k][k]
        a[k] = [x/t for x in a[k]]
        for j in range(n):
            if j != k:
                v = a[j][k]
                a[j] = [x-v*y for x, y in zip(a[j], a[k])]
    return [r[n:] for r in a]


def rowmul(v, a):
    return [sum(x*a[i][j] for i, x in enumerate(v)) for j in range(len(a))]


def basis(rows, n):
    """Integer row operations produce a triangular lattice basis."""
    assert all(Q(x).denominator == 1 for r in rows for x in r)
    a = [list(map(int, r)) for r in rows]
    for k in range(n):
        assert any(a[j][k] for j in range(k, len(a)))
        while True:
            j = min((j for j in range(k, len(a)) if a[j][k]), key=lambda j: abs(a[j][k]))
            a[j], a[k] = a[k], a[j]
            for j in range(k+1, len(a)):
                q = a[j][k]//a[k][k]
                a[j] = [x-q*y for x, y in zip(a[j], a[k])]
            if all(a[j][k] == 0 for j in range(k+1, len(a))):
                break
    assert all(not any(r) for r in a[n:])
    return a[:n]


class Field:
    def __init__(self, f, d):
        self.f, self.n = f, len(f)-1
        n = self.n
        assert factor_norms(f, 2) == [(2**n, 1)]
        self.eye = [[int(i == j) for j in range(n)] for i in range(n)]
        self.tracepowers = [Q(n)]
        for k in range(1, n):
            self.tracepowers.append(-sum(f[n-j]*self.tracepowers[k-j] for j in range(1, k))-k*f[n-k])
        self.gram = [[self.trace(self.multiply(a, b)) for a in self.eye] for b in self.eye]
        assert det(self.gram) == d

    def multiply(self, a, b):
        r = mul(a, b, self.f)
        return r + [Q(0)]*(self.n-len(r))

    def trace(self, a):
        return sum(x*y for x, y in zip(a, self.tracepowers))

    def norm(self, a):
        return det([self.multiply(a, e) for e in self.eye])

    def ideal_product(self, a, b):
        return basis([self.multiply(x, y) for x in a for y in b], self.n)

    def factors_of_element(self, a):
        norm = int(abs(self.norm(a)))
        assert norm > 0
        ps = [p for p in primes(norm) if norm % p == 0]
        out = []
        for p in ps:
            r = [x % p for x in self.f]
            # Only prime factors of norm <= norm(a) can divide (a).
            for degree in range(1, self.n+1):
                if p**degree > norm:
                    break
                for bits in product(range(p), repeat=degree):
                    g = list(bits)+[1]
                    if not any(divrem(r, g, p)[1]):
                        # Factors of smaller degree have already been removed.
                        count = 0
                        while not any(divrem(r, g, p)[1]):
                            r = divrem(r, g, p)[0]
                            count += 1
                        rows = [[p*x for x in e] for e in self.eye]
                        rows.extend(self.multiply(g, e) for e in self.eye)
                        P = basis(rows, self.n)
                        q = p**degree
                        assert abs(det(P)) == q
                        power = self.eye
                        v = 0
                        for _ in range(norm.bit_length()+1):
                            power = self.ideal_product(power, P)
                            member = all(x.denominator == 1 for x in rowmul(a, inverse(power)))
                            if not member:
                                break
                            v += 1
                        else:
                            raise AssertionError('valuation did not terminate')
                        if v:
                            out.append((q, v))
        assert prod(q**e for q, e in out) == norm
        return sorted(out)


def primes(bound):
    return [p for p in range(2, bound+1) if all(p % q for q in range(2, isqrt(p)+1))]


def val(f, x):
    y = Q(0)
    for a in reversed(f):
        y = y*x+a
    return y


def isolate(f):
    intervals = []
    M = 1+max(abs(x) for x in f[:-1])
    for j in range(-4*M, 4*M):
        a, b = Q(j, 4), Q(j+1, 4)
        if val(f, a)*val(f, b) < 0:
            for _ in range(75):
                m = (a+b)/2
                if val(f, m)*val(f, a) > 0:
                    a = m
                else:
                    b = m
            intervals.append((a, b))
    # Degree-many disjoint sign-changing intervals certify all roots.
    assert len(intervals) == len(f)-1
    return intervals


def intervalval(f, r):
    lo = hi = Q(0)
    for x in reversed(f):
        values = [lo*r[0], lo*r[1], hi*r[0], hi*r[1]]
        lo, hi = min(values)+x, max(values)+x
    return lo, hi


def positive(a, intervals):
    if not any(a):
        return False
    values = [intervalval(a, r) for r in intervals]
    assert all(lo > 0 or hi < 0 for lo, hi in values), 'refine root intervals'
    return all(lo > 0 for lo, hi in values)


def sqrtbound(x, upper=True):
    x = Q(x)
    den = 10**8
    k = isqrt((x.numerator*den**2)//x.denominator)
    return Q(k+1 if upper else k, den)


def verify(row):
    d = row['discriminant']
    f, A, weight, target = FIELDS[d]
    F = Field(f, d)
    n, q = F.n, 2**F.n
    roots = isolate(f)
    G = [[Q(x) for x in r] for r in row['lattice_certificate']['trace_gram_matrix']]
    assert G == F.gram
    inv = inverse(G)
    lengthbound = 16*n
    bounds = [isqrt(int(lengthbound*inv[i][i])) for i in range(n)]
    found = []
    tested = short = 0
    for c in product(*(range(-b, b+1) for b in bounds)):
        tested += 1
        length = sum(c[i]*G[i][j]*c[j] for i in range(n) for j in range(n))
        if length >= lengthbound:
            continue
        short += 1
        rest = [-x for x in c]
        rest[0] += 4
        if positive(c, roots) and positive(rest, roots):
            found.append(c)
    recorded = {tuple(z['coordinates']) for z in row['decompositions']}
    assert set(found) == recorded
    assert short == row['lattice_certificate']['number_of_short_vectors_including_zero']
    terms = []
    for z in row['decompositions']:
        x = list(z['coordinates'])
        y = [-v for v in x]
        y[0] += 4
        fx, fy = F.factors_of_element(x), F.factors_of_element(y)
        assert fx == sorted(map(tuple, z['ideal_x_factor_norms']))
        assert fy == sorted(map(tuple, z['ideal_y_factor_norms']))
        nx, ny = F.norm(x), F.norm(y)
        assert nx == z['norm_x'] and ny == z['norm_y']
        sigma = prod(sum(p**j for j in range(e+1)) for p, e in fx)
        tau = prod(e+1 for p, e in fy)
        assert sigma == z['sigma_one_x'] and tau == z['divisor_count_y']
        assert 0 < ny < q*q
        terms.append((sigma, tau, ny))
    alpha_certificate = None
    if d in (725, 5125):
        cert = row['alpha_upper_certificate']
        assert cert['cutoff'] == 199
        assert [e['p'] for e in cert['Euler_factors']] == primes(199)
        P = Q(1)
        for e in cert['Euler_factors']:
            fact = factor_norms(f, e['p'])
            assert fact == sorted((v['norm'], v['ramification_index']) for v in e['primes'])
            local = prod(Q(p*p, p*p-1) for p, exponent in fact)
            assert local == Q(e['factor'])
            P *= local
        assert P == Q(cert['Euler_product'])
        alpha_upper = (4*Q(22, 7)**2)**n/(d*sqrtbound(d, False)*P)
        assert alpha_upper < A
        alpha_certificate = {'Euler_product': str(P), 'alpha_upper': str(alpha_upper), 'A': str(A)}
    v = weight-1
    rhs = A*(sum(Q(sigma*tau)*(sqrtbound(ny)/q)**v for sigma, tau, ny in terms)
             +4*sqrtbound(Q(1, q**v))+A/q**v)
    assert rhs < target < q*q-1
    return {'discriminant': d, 'degree': n, 'coordinate_points_checked': tested,
            'short_vectors': short, 'positive_decompositions': len(found),
            'source_weight_excluded_from': weight, 'A': str(A),
            'rational_rhs_upper': str(rhs), 'integer_upper': target,
            'left_side': q*q-1, 'Euler_alpha_certificate': alpha_certificate,
            'larger_weights_excluded_by': 'all exact norm bases lie strictly between zero and one'}


def main():
    if not __debug__:
        raise SystemExit('Assertions must be enabled: do not use -O/-OO.')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parent/'rerun')
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    path = args.data_dir/'rerun'/'e2_dyadic_output.json'
    if not path.is_file():
        path = args.data_dir/'e2_dyadic_output.json'
    data = json.loads(path.read_text())
    rows = [r for r in data['fields'] if r.get('discriminant') in FIELDS]
    assert len(rows) == 6 and {r['discriminant'] for r in rows} == set(FIELDS)
    result = {'status': 'passed', 'arithmetic': 'integer/rational, Python standard library',
              'input_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'external_dependencies': ['stated field discriminants', 'four exact alpha values',
                                        'coefficient inequality and Ramanujan bound'],
              'fields': [verify(r) for r in rows]}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: six exhaustive dyadic certificates and weight cutoffs')


if __name__ == '__main__':
    main()
