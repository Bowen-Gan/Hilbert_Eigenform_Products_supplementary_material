#!/usr/bin/env python3
"""Run twelve original Sage tasks followed by six additional exact checks.

The default invocation records both suites and reports complete success only
when both pass in this invocation. Endpoint-only mode retains its three tasks.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import uuid

ROOT = Path(__file__).resolve().parent


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    return json.loads(path.read_text())


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2)+'\n')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes():
    return {p.name: digest(p) for p in sorted(ROOT.iterdir())
            if p.is_file() and p.suffix in {'.py', '.sage'}}


def refresh_checksums():
    files = []
    for directory, subdirs, names in os.walk(ROOT):
        subdirs[:] = sorted(d for d in subdirs if d not in {'.git', '__pycache__'})
        for name in names:
            path = Path(directory)/name
            if path == ROOT/'SHA256SUMS.txt' or path.suffix in {'.pyc', '.pyo'}:
                continue
            files.append(path)
    lines = [digest(p)+'  '+p.relative_to(ROOT).as_posix()
             for p in sorted(files, key=lambda p: p.relative_to(ROOT).as_posix())]
    (ROOT/'SHA256SUMS.txt').write_text('\n'.join(lines)+'\n')


def additional_checks():
    result = subprocess.run(
        [sys.executable, str(ROOT/'run_verification.py'),
         '--data-dir', str(ROOT/'rerun'),
         '--output-dir', str(ROOT/'verification_case11')],
        cwd=ROOT, env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
    if result.returncode:
        return result.returncode
    data = read(ROOT/'verification_case11'/'status.json')
    assert data['completed'] and len(data['tasks']) == 6
    assert all(row['exit_code'] == 0 and row['accepted'] for row in data['tasks'])
    assert all(data['cross_checks'].values())
    assert all(digest(ROOT/name) == checksum
               for name, checksum in data['source_hashes'].items())
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--endpoints-only', action='store_true')
    modes.add_argument('--additional-only', action='store_true')
    args = parser.parse_args()
    if sys.flags.optimize:
        raise RuntimeError('Assertions must remain enabled; do not use -O/-OO.')
    sys.dont_write_bytecode = True
    if args.additional_only:
        result = additional_checks()
        if result == 0:
            refresh_checksums()
        return result

    run_id = str(uuid.uuid4())
    full_path = ROOT/'verification'/'full_run_status.json'
    full = {
        'schema': 'complete-arithmetic-run-v1', 'run_id': run_id,
        'scope': 'twelve original Sage tasks plus six additional exact checks',
        'started_utc': now(), 'completed': False, 'stage': 'sage_environment',
        'python_version': sys.version, 'platform': platform.platform(),
        'source_hashes': source_hashes(), 'suites': [],
        'theoretical_inputs': 'The manuscript and cited mathematical theorems remain external.',
    }
    if not args.endpoints_only:
        write(full_path, full)
    try:
        import sage.all
        from sage.version import version as sage_version
    except ImportError:
        if not args.endpoints_only:
            full.update(stage='failed', error='SageMath is unavailable in this Python environment.',
                        finished_utc=now())
            write(full_path, full)
        print('SageMath is unavailable. Activate your Sage Python environment and rerun.',
              file=sys.stderr)
        return 1

    out = ROOT/'rerun'
    out.mkdir(exist_ok=True)
    def output(name):
        return ['--output', str(out/name)]
    tasks = [] if args.endpoints_only else [
        ('e3_cusp_check.sage', []),
        ('e2_square_branch_check.sage', []),
        ('e2_high_degree_check.sage', []),
        ('e2_candidate_enumeration.sage', output('e2_candidates.json')),
        ('e2_allclass_balanced.sage', output('d3969_certificate.json')),
        ('e2_indecomposable_search.sage',
         ['--input', str(out/'e2_candidates.json'), *output('e2_lattices.json')]),
        ('e2_dyadic_bound.sage', output('e2_dyadic_output.json')),
        ('rq_certificate.py', output('rq_certificate.json')),
        ('rq_e2_exact.sage', []),
    ]
    tasks += [
        ('cubic_weight_two_genus.sage', output('cubic_weight_two_genus.json')),
        ('quartic_725_weight_two.sage', output('quartic_725_weight_two.json')),
        ('quartic_1125_sufficient.sage', output('quartic_1125_sufficient.json')),
    ]
    scope = 'endpoints' if args.endpoints_only else 'all_arithmetic'
    status = {
        'scope': scope, 'run_id': run_id, 'started_utc': now(),
        'completed': False, 'tasks': [], 'sage_version': sage_version,
        'python_version': sys.version, 'platform': platform.platform(),
        'source_hashes': source_hashes(), 'magma_required': False,
        'theoretical_inputs': 'Read notes/; interpretations use the manuscript and cited theorems.',
    }
    status_path = out/(scope+'_status.json')
    write(status_path, status)
    if not args.endpoints_only:
        full.update(stage='original_arithmetic', sage_version=sage_version)
        write(full_path, full)
    try:
        for name, options in tasks:
            print('RUN '+name, flush=True)
            log = out/(Path(name).stem+'.txt')
            command = [sys.executable, str(ROOT/name), *options]
            with log.open('w') as handle:
                result = subprocess.run(command, cwd=ROOT, stdout=handle,
                                        stderr=subprocess.STDOUT,
                                        env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
            accepted = result.returncode == 0
            if name == 'e2_indecomposable_search.sage' and result.returncode == 2:
                data = read(out/'e2_lattices.json')
                accepted = (data['search_complete'] and data['all_required_degrees']
                            and sorted(r['discriminant'] for r in data['unresolved_fields'])
                            == [169, 361, 725])
            status['tasks'].append({'file': name, 'command': command,
                                    'exit_code': result.returncode, 'accepted': bool(accepted),
                                    'log': log.name, 'log_sha256': digest(log)})
            status['updated_utc'] = now()
            write(status_path, status)
            if not accepted:
                raise RuntimeError('Original arithmetic task failed: '+name+'; inspect '+str(log))

        cubic = read(out/'cubic_weight_two_genus.json')
        assert {r['discriminant'] for r in cubic['rows']} == {169, 361}
        assert all(r['full_level_parallel_weight_two_cusp_dimension'] == 0 for r in cubic['rows'])
        assert read(out/'quartic_725_weight_two.json')['full_level_parallel_weight_two_cusp_dimension'] == 0
        large = read(out/'quartic_1125_sufficient.json')
        assert large['conclusions']['dim_S2_trivial_character'] == 0
        assert large['conclusions']['dim_S5_totally_odd_character_at_least'] >= 5
        if not args.endpoints_only:
            enumeration = read(out/'e2_candidates.json')
            assert enumeration['covers_all_required_degrees'] and len(enumeration['fields']) == 771
            lattices = read(out/'e2_lattices.json')
            assert sum(r['status'] == 'excluded_by_lattice_certificate' for r in lattices['fields']) == 15
        assert source_hashes() == status['source_hashes']
        status.update(completed=True, finished_utc=now(),
                      new_arithmetic_endpoints_closed=[169, 361, 725, 1125])
        status['output_hashes'] = {p.name: digest(p) for p in sorted(out.glob('*.json'))
                                   if p != status_path and not p.name.endswith('_status.json')}
        write(status_path, status)
        if args.endpoints_only:
            refresh_checksums()
            print('FREE_ENDPOINTS_PASSED')
            return 0

        full['suites'].append({'name': 'original_arithmetic', 'task_count': 12,
                               'completed': True, 'status': 'rerun/all_arithmetic_status.json',
                               'status_sha256': digest(status_path)})
        full['stage'] = 'additional_case11'
        write(full_path, full)
        additional_started = now()
        result = additional_checks()
        if result:
            raise RuntimeError('Additional Case 11 checks failed; inspect verification_case11/.')
        extra_path = ROOT/'verification_case11'/'status.json'
        extra = read(extra_path)
        assert extra['started_utc'] >= additional_started
        assert extra['input_hashes'] == {name: digest(out/name) for name in
                                       ['e2_lattices.json', 'd3969_certificate.json', 'e2_dyadic_output.json']}
        assert source_hashes() == full['source_hashes']
        full['suites'].append({'name': 'additional_case11', 'task_count': 6,
                               'completed': True, 'status': 'verification_case11/status.json',
                               'status_sha256': digest(extra_path)})
        full.update(completed=True, stage='completed', finished_utc=now(),
                    success_marker='FULL_CERTIFICATES_PASSED')
        write(full_path, full)
        refresh_checksums()
        print('FREE_ARITHMETIC_SUITE_PASSED')
        print('FULL_CERTIFICATES_PASSED')
        return 0
    except Exception as exc:
        if not args.endpoints_only:
            full.pop('success_marker', None)
            full.update(completed=False, stage='failed', error=str(exc), finished_utc=now())
            write(full_path, full)
        print('FAILED: '+str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
