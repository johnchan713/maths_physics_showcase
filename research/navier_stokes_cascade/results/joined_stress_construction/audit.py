#!/usr/bin/env python3
"""Reproduce the construction proposal and the review that withholds promotion."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import mpmath as mp
from construction import (HERE, RESULTS, CANDIDATE, OPEN_OBLIGATIONS, geometry,
                          derivative_constants, c1_counterexamples,
                          finite_integer_comparison, promotion_review, lower, upper)
from review import curvature_identity, loop_derivative_controls, rounded_endpoint_control

STATUS = 'joined-stress-proposal-reviewed-derivative-budget-open'


def source_hashes():
    paths = [HERE/n for n in ('README.md', 'protocol.json', 'requirements.txt',
                              'construction.py', 'review.py', 'audit.py', 'test_audit.py')]
    for directory, names in {
        'axis_core_attachment': ('README.md', 'bounds.py', 'evidence.json'),
        'axis_matching_audit': ('README.md', 'moments.py', 'evidence.json'),
        'stress_realization_audit': ('README.md', 'loop.py', 'bounds.py', 'repair.py', 'evidence.json'),
        'intermediate_decay_audit': ('evidence.json',),
        'heat_exterior_audit': ('README.md', 'evidence.json'),
        'outer_pressure_pilot': ('schedule.py',),
    }.items():
        paths += [RESULTS/directory/n for n in names]
    return {str(path.relative_to(RESULTS.parent)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def inherited_hashes():
    results = {}
    for directory in ('axis_core_attachment', 'axis_matching_audit', 'stress_realization_audit'):
        data = json.loads((RESULTS/directory/'evidence.json').read_text())
        hashes = data.get('source_hashes', data.get('provenance'))
        if not isinstance(hashes, dict) or not hashes:
            raise ValueError('Missing inherited hashes: '+directory)
        results[directory] = all(hashlib.sha256((RESULTS.parent/path).read_bytes()).hexdigest() == sha
                                 for path, sha in hashes.items())
    return results


def encode(value):
    if hasattr(value, '_mpi_'):
        return dict(lower=mp.nstr(lower(value), 120), upper=mp.nstr(upper(value), 120))
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 90)
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    return value


def generate():
    mp.mp.dps = 80
    protocol = json.loads((HERE/'protocol.json').read_text())
    if protocol['status'] != STATUS or protocol['interval_digits'] != [80, 110] \
            or protocol['actual_joined_profile_frequency_selected'] is not False:
        raise ValueError('Changed protocol or unsupported scientific promotion')
    geometries = [geometry(d) for d in protocol['interval_digits']]
    derivatives = [derivative_constants(d) for d in protocol['interval_digits']]
    curvature, loop = curvature_identity(), loop_derivative_controls()
    endpoints, integers, inherited = rounded_endpoint_control(), finite_integer_comparison(), inherited_hashes()
    controls = c1_counterexamples()
    eps = Fraction(1, 10**16)
    gates = dict(
        all_actual_geometry_checks=all(all(g['checks'].values()) for g in geometries),
        all_universal_derivative_constants=all(all(g['checks'].values()) for g in derivatives),
        exact_polynomial_cap_comparisons=all(integers.values()),
        second_derivative_identity_and_missing_term_controls=all(curvature['checks'].values()),
        signed_variance_derivative_controls=loop['check'],
        narrow_endpoint_rounding_failure_retained=endpoints['rounded_offsets_collapse'],
        C1_smallness_does_not_bound_C2=all(r['C1_upper'] <= eps for r in controls)
            and controls[-1]['second_derivative_sup'] > 1,
        frozen_parent_sources_unchanged=all(inherited.values()),
        unproved_actual_frequency_rejected=promotion_review()['accepted'] is False,
    )
    return dict(status=STATUS, protocol=protocol, source_hashes=source_hashes(),
                geometry=geometries, derivative_transfer=derivatives,
                implicit_second_derivative_review=curvature,
                signed_variance_derivative_review=loop,
                endpoint_rounding_review=endpoints, C1_counterexamples=controls,
                cap_comparisons=integers, proposed_construction=CANDIDATE,
                open_obligations=OPEN_OBLIGATIONS, promotion=promotion_review(),
                inherited_provenance=inherited, gates=gates)


def validate(record):
    if record.get('status') != STATUS or record.get('promotion') != promotion_review():
        raise ValueError('The current review does not certify an actual-profile frequency or a blow-up')
    if record.get('open_obligations') != OPEN_OBLIGATIONS:
        raise ValueError('An unresolved estimate was removed or changed')
    if record.get('proposed_construction') != CANDIDATE:
        raise ValueError('The construction proposal was changed without a new review')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A construction/review check failed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-record', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    record = encode(generate())
    validate(record)
    if args.verify_record and json.loads(args.verify_record.read_text()) != record:
        raise SystemExit('Frozen evidence differs; no record was written')
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    print(STATUS)
    print(f"{len(record['gates'])} construction/review gates passed; frequency promotion rejected.")
    print('Actual derivative transfer, pressure-error constants, full PDE and blow-up remain unverified.')


if __name__ == '__main__':
    main()
