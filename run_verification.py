#!/usr/bin/env python3
"""Run ten independent exact-arithmetic verifiers and archive this run.

No Sage installation is needed. This is an additional verification suite,
not a substitute for the original Sage field-enumeration and endpoint programs.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import sys
import uuid


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if not __debug__:
        raise SystemExit('Do not use optimized Python; assertions are required.')
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-dir', type=Path, default=root/'rerun')
    ap.add_argument('--output-dir', type=Path, default=root/'verification_case11')
    ap.add_argument('--parent-run-id', help=argparse.SUPPRESS)
    args = ap.parse_args()
    inputs, out = args.data_dir.resolve(), args.output_dir.resolve()
    complete_path = root/'verification'/'full_run_status.json'
    if args.parent_run_id:
        complete = json.loads(complete_path.read_text())
        if complete['run_id'] != args.parent_run_id or complete['completed']:
            raise RuntimeError('Independent suite does not belong to the active complete run.')
    else:
        partial = {
            'schema': 'complete-arithmetic-run-v2', 'run_id': str(uuid.uuid4()),
            'scope': 'independent checks only; not a complete reproduction',
            'completed': False, 'stage': 'independent_only',
            'started_utc': datetime.now(timezone.utc).isoformat(),
        }
        complete_path.parent.mkdir(parents=True, exist_ok=True)
        complete_path.write_text(json.dumps(partial, indent=2)+'\n')
        (root/'STATUS.json').write_text(json.dumps(partial, indent=2)+'\n')
    out.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent
    selected_inputs = {
        'e2_lattices.json': inputs/'e2_lattices.json',
        'd3969_certificate.json': inputs/'d3969_certificate.json',
        'e2_dyadic_output.json': (inputs/'rerun'/'e2_dyadic_output.json'
                                 if (inputs/'rerun'/'e2_dyadic_output.json').is_file()
                                 else inputs/'e2_dyadic_output.json'),
        'e2_candidates.json': inputs/'e2_candidates.json',
        'quartic_1125_sufficient.json': inputs/'quartic_1125_sufficient.json',
    }
    tasks = [
        ('verify_lattice_witnesses.py', ['--data-dir', str(inputs)], 'lattice_verification.json'),
        ('verify_dyadic_bounds.py', ['--data-dir', str(inputs)], 'dyadic_verification.json'),
        ('verify_special_values.py', [], 'special_values_verification.json'),
        ('verify_d21_local.py', [], 'd21_d40_local_result.json'),
        ('rq_certificate.py', [], 'quadratic_reduction_verification.json'),
        ('verify_rankin_local_factor.py', [], 'rankin_local_factor.json'),
        ('verify_field_enumeration.py', ['--data-dir', str(inputs)], 'field_enumeration_verification.json'),
        ('verify_quaternion_central_character.py', ['--data-dir', str(inputs)], 'quaternion_central_character.json'),
        ('verify_d12_weight_one.py', [], 'd12_weight_one_verification.json'),
        ('verify_eisenstein_bounds.py', [], 'eisenstein_bounds_verification.json'),
    ]
    record = {
        'schema': 'independent-arithmetic-run-v2',
        'parent_run_id': args.parent_run_id,
        'scope': 'additional arithmetic checks; analytic theorems and complete field enumeration are external',
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'python_version': sys.version, 'platform': platform.platform(),
        'completed': False, 'tasks': [],
        'source_hashes': {p.name: digest(p) for p in sorted(scripts.glob('*.py'))},
        'input_paths': {name: str(p) for name, p in selected_inputs.items()},
        'input_hashes': {},
    }
    status = out/'status.json'
    status.write_text(json.dumps(record, indent=2)+'\n')
    missing = [str(p) for p in selected_inputs.values() if not p.is_file()]
    if missing:
        record['error'] = 'Required input files are missing: '+', '.join(missing)
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        status.write_text(json.dumps(record, indent=2)+'\n')
        raise SystemExit(record['error'])
    record['input_hashes'] = {name: digest(p) for name, p in selected_inputs.items()}
    status.write_text(json.dumps(record, indent=2)+'\n')
    for name, extra, result_name in tasks:
        print('RUN '+name, flush=True)
        (out/result_name).unlink(missing_ok=True)
        command = [sys.executable, str(scripts/name), *extra, '--output', str(out/result_name)]
        run = subprocess.run(command, capture_output=True, text=True,
                             cwd=root, env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
        log = out/(Path(name).stem+'.log')
        log.write_text(run.stdout+run.stderr)
        row = {'script': name, 'command': command, 'exit_code': run.returncode,
               'accepted': run.returncode == 0, 'log': log.name, 'result': result_name}
        record['tasks'].append(row)
        record['updated_utc'] = datetime.now(timezone.utc).isoformat()
        status.write_text(json.dumps(record, indent=2)+'\n')
        if run.returncode:
            print('FAILED '+name+'; inspect '+str(log), file=sys.stderr)
            return 1
        result = out/result_name
        if not result.is_file():
            row['accepted'] = False
            row['error'] = 'expected result missing'
            status.write_text(json.dumps(record, indent=2)+'\n')
            return 1
        json.loads(result.read_text())
        row['result_sha256'] = digest(result)
    # Connect the independently computed constants to the bound inputs.
    special = json.loads((out/'special_values_verification.json').read_text())
    dyadic = json.loads((out/'dyadic_verification.json').read_text())
    constants = {r['discriminant']: r['alpha'] for r in special['fields']}
    from fractions import Fraction
    for r in dyadic['fields']:
        if r['discriminant'] in constants:
            assert Fraction(r['A']) == abs(Fraction(constants[r['discriminant']]))
    quadratic = json.loads((out/'quadratic_reduction_verification.json').read_text())
    assert quadratic['reduction']['survivors'] == [12,21,24,28,69,77]
    assert next(r for r in quadratic['reduction']['rows'] if r['D'] == 21)['alpha'] == ['12']
    rankin = json.loads((out/'rankin_local_factor.json').read_text())
    assert rankin['status'] == 'PASS' and rankin['difference'] == []
    record['cross_checks'] = {
        'special_values_match_four_dyadic_inputs': True,
        'quadratic_reduction_has_exactly_six_survivors': True,
        'D21_reciprocal_constant_is_12': True,
        'Rankin_local_identity_is_exact': True,
    }
    fields = json.loads((out/'field_enumeration_verification.json').read_text())
    quaternion = json.loads((out/'quaternion_central_character.json').read_text())
    weight_one = json.loads((out/'d12_weight_one_verification.json').read_text())
    eisenstein = json.loads((out/'eisenstein_bounds_verification.json').read_text())
    assert fields['status'] == 'PASS'
    assert quaternion['result'] == 'PASS'
    assert weight_one['status'] == 'PASS'
    assert eisenstein['status'] == 'passed' and eisenstein['completed'] and eisenstein['passed']
    record['cross_checks'].update({
        'complete_number_field_inputs_match_recomputed_enumeration': True,
        'quaternion_invariants_have_required_central_character': True,
        'weight_one_dimensions_and_normalization_are_exact': True,
        'finite_eisenstein_product_bounds_are_certified': True,
    })
    record['completed'] = True
    record['finished_utc'] = datetime.now(timezone.utc).isoformat()
    record['success_marker'] = 'CASE11_ADDITIONAL_ARITHMETIC_PASSED'
    status.write_text(json.dumps(record, indent=2)+'\n')
    print(record['success_marker'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
