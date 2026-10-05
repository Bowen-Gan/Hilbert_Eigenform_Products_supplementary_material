#!/usr/bin/env python3
"""Independent exact check of the 772-field arithmetic-screening manifest.

Python standard library only.  Check every bundled source row and hash, the
complete cutoff coverage, source-to-output identities, trace discriminants,
integer-order bases and ring closure, power-basis indices, ideal norms,
dyadic ramification via Frobenius on O/2O, exact Euler products, and rational
intervals obtained independently from Machin's formula for pi.

External inputs remain explicit: completeness and field discriminants of
Bordeaux's published tables; ordinary/narrow class numbers and primality of
non-dyadic local factors, recomputed by the separate proof-enabled Sage run;
the manuscript's modular-form exclusions. This program verifies their recorded
arithmetic use; it does not claim to reprove field enumeration or class groups.

Run from the supplement root:
    python verify_field_enumeration.py --data-dir rerun --output verification_case11/field_enumeration.json
"""
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt, prod
from collections import Counter
import argparse
import hashlib
import json
import sys
from number_field_table_inputs import ENUMERATION_BOUNDS, EXPECTED_COUNTS, table_rows, load_manifest

if not __debug__:
    raise RuntimeError("verification assertions require Python without -O/-OO")


def determinant(matrix):
    a = [list(row) for row in matrix]
    value = Q(1)
    for k in range(len(a)):
        if not any(a[j][k] for j in range(k, len(a))):
            return Q(0)
        j = next(j for j in range(k, len(a)) if a[j][k])
        if j != k:
            a[k], a[j] = a[j], a[k]
            value = -value
        pivot = a[k][k]
        value *= pivot
        for j in range(k + 1, len(a)):
            scale = a[j][k]/pivot
            for t in range(k + 1, len(a)):
                a[j][t] -= scale*a[k][t]
    return value


def inverse(matrix):
    n = len(matrix)
    a = [list(row) + [Q(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for k in range(n):
        j = next(j for j in range(k, n) if a[j][k])
        a[k], a[j] = a[j], a[k]
        pivot = a[k][k]
        a[k] = [v/pivot for v in a[k]]
        for j in range(n):
            if j != k:
                scale = a[j][k]
                a[j] = [v-scale*w for v, w in zip(a[j], a[k])]
    return [row[n:] for row in a]


def row_times(v, matrix):
    return [sum(v[i]*matrix[i][j] for i in range(len(v))) for j in range(len(v))]


def polynomial_product(a, b, f):
    n = len(f)-1
    c = [Q(0)]*(2*n-1)
    for i in range(n):
        for j in range(n):
            c[i+j] += a[i]*b[j]
    for k in range(2*n-2, n-1, -1):
        for j in range(n):
            c[k-n+j] -= c[k]*f[j]
    return c[:n]


def basis_matrix(rows, n):
    assert len(rows) == n
    return [[Q(v) for v in row] + [Q(0)]*(n-len(row)) for row in rows]


def rank_mod_two(matrix):
    a = [[int(v) % 2 for v in row] for row in matrix]
    rank = 0
    for k in range(len(a[0])):
        candidates = [i for i in range(rank, len(a)) if a[i][k]]
        if not candidates:
            continue
        i = candidates[0]
        a[rank], a[i] = a[i], a[rank]
        for i in range(len(a)):
            if i != rank and a[i][k]:
                a[i] = [x ^ y for x, y in zip(a[i], a[rank])]
        rank += 1
    return rank


def atan_bounds(q, terms=50):
    s = sum((Q(1) if j % 2 == 0 else Q(-1)) /
            ((2*j+1)*q**(2*j+1)) for j in range(terms))
    next_term = (Q(1) if terms % 2 == 0 else Q(-1)) / ((2*terms+1)*q**(2*terms+1))
    return min(s, s+next_term), max(s, s+next_term)


def pi_bounds():
    a, b = atan_bounds(5), atan_bounds(239)
    lo, hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    denominator = 10**60
    lo = Q((lo*denominator).__floor__(), denominator)
    hi = Q((hi*denominator).__ceil__(), denominator)
    assert Q(3141592653589793238462643383279, 10**30) < lo < hi
    return lo, hi


def primes_through(bound):
    return [p for p in range(2, bound+1) if all(p % d for d in range(2, isqrt(p)+1))]


def positive_interval_power(interval, exponent):
    return interval[0]**exponent, interval[1]**exponent


def certify_order(row):
    n, f = row['degree'], list(map(Q, row['polynomial']))
    assert len(f) == n+1 and f[-1] == 1
    traces = [Q(n)]
    for k in range(1, 2*n-1):
        value = -sum(f[n-j]*traces[k-j] for j in range(1, min(k, n)+1))
        if k <= n:
            value += (n-k)*f[n-k]  # replace f[n-k]*Tr(1) by k*f[n-k]
        traces.append(value)
    gram = [[traces[i+j] for j in range(n)] for i in range(n)]
    power_discriminant = determinant(gram)
    assert power_discriminant == row['power_basis_discriminant']
    index = row['power_basis_index']
    assert index >= 1 and power_discriminant == row['discriminant']*index**2
    order = basis_matrix(row['maximal_order_basis'], n)
    inv = inverse(order)
    assert abs(determinant(order)) == Q(1, index)
    assert determinant(order)**2*power_discriminant == row['discriminant']
    assert all(v.denominator == 1 for v in row_times([Q(1)]+[Q(0)]*(n-1), inv))
    squares = []
    for i in range(n):
        for j in range(n):
            coordinates = row_times(polynomial_product(order[i], order[j], f), inv)
            assert all(v.denominator == 1 for v in coordinates)
            if i == j:
                squares.append(coordinates)
    frobenius_rank = rank_mod_two(squares)
    dyadic = row['dyadic_factorization']
    assert sum(v['ramification_index']*v['residue_degree'] for v in dyadic) == n
    assert (frobenius_rank < n) == any(v['ramification_index'] > 1 for v in dyadic)
    check_local_bases(row, order, inv, f, 2, dyadic)
    for local in row.get('small_prime_factorizations', []):
        check_local_bases(row, order, inv, f, local['rational_prime'], local['factors'])
    return frobenius_rank


def check_local_bases(row, order, order_inverse, f, p, factors):
    n = row['degree']
    assert sum(v['ramification_index']*v['residue_degree'] for v in factors) == n
    for factor in factors:
        ideal = basis_matrix(factor['basis'], n)
        coordinates = [row_times(v, order_inverse) for v in ideal]
        assert all(v.denominator == 1 for rowv in coordinates for v in rowv)
        assert abs(determinant(coordinates)) == factor['norm'] == p**factor['residue_degree']
        inv = inverse(ideal)
        for a in order:
            assert all(v.denominator == 1 for v in row_times([p*c for c in a], inv))
            for b in ideal:
                assert all(v.denominator == 1 for v in row_times(polynomial_product(a, b, f), inv))


def check_euler(row, pi):
    data = row['alpha_interval']
    cutoff, n, D = data['cutoff'], row['degree'], row['discriminant']
    assert cutoff == 199
    assert [p for p, degrees in data['factor_degrees']] == primes_through(cutoff)
    P = Q(1)
    for p, degrees in data['factor_degrees']:
        assert all(type(d) is int and d >= 1 for d in degrees)
        assert sum(degrees) <= n
        if D % p != 0:
            assert sum(degrees) == n
        for d in degrees:
            P /= 1-Q(1, p**(2*d))
    assert P == Q(data['euler_product'])
    tail = Q(cutoff+1, cutoff)**n
    assert tail == Q(data['tail_bound'])
    scale = 10**60
    root = isqrt(D*scale**2)
    sqrt_lo, sqrt_hi = Q(root, scale), Q(root+1, scale)
    lower = (4*pi[0]**2)**n / (D*sqrt_hi*P*tail)
    upper = (4*pi[1]**2)**n / (D*sqrt_lo*P)
    integer_min, integer_max = lower.__ceil__(), upper.__floor__()
    assert integer_min == data['possible_integer_min']
    assert integer_max == data['possible_integer_max']
    assert data['excludes_integrality'] == (integer_min > integer_max)
    return integer_min > integer_max


def check_status(row, pi):
    n, h, D = row['degree'], row['ordinary_class_number'], row['discriminant']
    assert row['class_number_proof'] is True
    if h % 2 == 0:
        expected = 'even_class_number_argument'
    elif n*h > 14:
        expected = 'general_degree_class_number_bound'
    else:
        lower = (4*pi[0]**2)**(2*n)
        upper = (4*pi[1]**2)**(2*n)
        assert not lower <= D**3 <= upper
        if h == 1 and D**3 > upper:
            expected = 'trivial_character_discriminant_bound'
        else:
            factors = row['dyadic_factorization']
            inert = len(factors) == 1 and factors[0]['ramification_index'] == 1 and factors[0]['residue_degree'] == n
            if not inert:
                expected = 'noninert_two_argument'
            else:
                tested = row['small_prime_factorizations']
                bad = [item for local in tested for item in local['factors'] if item['norm']**2 <= 2**n]
                if bad:
                    expected = 'small_prime_argument'
                    witness = row['exclusion_prime']
                    assert witness['p'] in primes_through(isqrt(2**n)) and witness['p'] > 2
                    assert any(local['rational_prime'] == witness['p'] and
                               any(v['norm'] == witness['norm'] and v['basis'] == witness['basis'] for v in local['factors'])
                               for local in tested)
                else:
                    assert [local['rational_prime'] for local in tested] == primes_through(isqrt(2**n))[1:]
                    nonintegral = check_euler(row, pi) if h == 1 else False
                    expected = 'reciprocal_constant_not_integral' if nonintegral else 'requires_lattice_or_exception_certificate'
    assert row['status'] == expected, (D, row['status'], expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parent/'rerun')
    parser.add_argument('--output', type=Path, default=Path('verification_case11/field_enumeration.json'))
    args = parser.parse_args()
    path = args.data_dir/'e2_candidates.json'
    raw = path.read_bytes()
    data = json.loads(raw)
    assert data['enumeration_complete'] and data['covers_all_required_degrees']
    assert set(data['requested_degrees']) == set(ENUMERATION_BOUNDS)
    assert data['endpoint_checks'] == {str(k): v for k, v in ENUMERATION_BOUNDS.items()}
    assert data['input_provenance'] == load_manifest()
    pi = pi_bounds()
    source = []
    assert len(data['requested_degrees']) == len(ENUMERATION_BOUNDS)
    for n in data['requested_degrees']:
        cutoff = ENUMERATION_BOUNDS[n]
        base = (2*pi[0]**4/3, 2*pi[1]**4/3) if n in (3, 4) else (4*pi[0]**2, 4*pi[1]**2)
        assert (cutoff-1)**3 < base[0]**(2*n) < base[1]**(2*n) < cutoff**3
        source.extend(table_rows(n, cutoff))
    # The septic source completeness/local boundary is regenerated separately.
    assert len(table_rows(7)) == 4
    rows = data['fields']
    assert len(rows) == len(source) == 772
    ranks = {}
    for row, original in zip(rows, source):
        for key in ('degree', 'discriminant', 'polynomial'):
            assert row[key] == original[key]
        assert row['ordinary_class_number'] == original['table_ordinary_class_number']
        assert row['source_table'] == {k: original[k] for k in ('table_file', 'table_line', 'table_ordinary_class_number', 'table_class_group_invariants', 'table_polynomial_highest_degree_first')}
        ranks[(row['degree'], row['discriminant'], tuple(row['polynomial']))] = certify_order(row)
        check_status(row, pi)
    counts = Counter(str(row['degree']) for row in rows)
    statuses = Counter(row['status'] for row in rows)
    assert counts == data['counts_by_degree'] == {str(n): EXPECTED_COUNTS[n] for n in ENUMERATION_BOUNDS}
    assert statuses == data['counts_by_status']
    survivors = [row for row in rows if row['status'] == 'requires_lattice_or_exception_certificate']
    assert len(survivors) == 21
    assert Counter(row['degree'] for row in survivors) == {3: 10, 4: 8, 5: 2, 6: 1}
    omitted = [row for row in rows if row['discriminant'] == 65808]
    assert len(omitted) == 1 and omitted[0]['power_basis_index'] == 8
    key = (4, 65808, tuple(omitted[0]['polynomial']))
    assert ranks[key] == 2
    assert sorted((v['ramification_index'], v['residue_degree'], v['norm'])
                  for v in omitted[0]['dyadic_factorization']) == [(2, 1, 2), (2, 1, 2)]
    payload = {'schema': 'case11-independent-field-screening-v1', 'result': 'PASS',
               'status': 'PASS',
               'candidate_manifest_sha256': hashlib.sha256(raw).hexdigest(),
               'field_count': len(rows), 'counts_by_degree': dict(counts),
               'counts_by_status': dict(statuses), 'remaining_field_count': len(survivors),
               'D65808_power_basis_index': 8,
               'D65808_frobenius_rank_mod_two': ranks[key],
               'D65808_ramification_independently_checked': True,
               'external_inputs': __doc__.split('External inputs remain explicit:')[1].split('Run from')[0].strip()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2)+'\n')
    print('ALL_PASS: fields=772 degrees=143,552,37,40 remaining=21 D65808_index=8 ramified=True')


if __name__ == '__main__':
    main()
