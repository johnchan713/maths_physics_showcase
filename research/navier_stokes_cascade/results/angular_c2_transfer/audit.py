#!/usr/bin/env python3
"""Reproduce the actual C2 transfer argument's constants, controls, and scope."""
import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import mpmath as mp
from bounds import (HERE, RESULTS, previous, core_transfer, normalization_constants,
                    matching_transfer, lower, upper)
from review import (normalization_review, cancellation_review, exponential_review,
                    axis_integrability_review)

STATUS = 'actual-C2-moment-transfer-and-matching-bounded'
FLAGS = dict(actual_target_C2_transfer_bounded=True,
             actual_matching_coefficients_C2_bounded=True,
             same_exact_joining_correction=True,
             relies_on_preceding_analytic_core_argument=True,
             actual_compact_jet_envelope_verified=False,
             actual_joined_profile_frequency_selected=False,
             full_admissible_stress_realized=False,
             full_PDE_corrections_verified=False,
             smooth_force_verified=False, blowup_verified=False,
             formal_proof_assistant_certificate=False,
             independent_peer_review_completed=False)
RESOLVED = ['actual_target_C2_transfer']
REMAINING = {k: v for k, v in previous.OPEN_OBLIGATIONS.items() if k not in RESOLVED}


def validate_protocol(protocol):
    expected = dict(status=STATUS, parent_commit='86f5b77183ed9f93f8faa3fd4bee548596bf5b53',
                    source_pdf_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                    source_pages=[145, 146, 147, 148, 150, 154, 155, 156],
                    scale='Md64', interval_digits=[80, 110],
                    angular_norm='sup(abs(f)) + sup(abs(f_eta)) + sup(abs(f_etaeta))/2',
                    entry_epsilon='1e-16', target_C2_bound_over_epsilon='0.286',
                    actual_joined_profile_frequency_selected=False)
    if protocol != expected:
        raise ValueError('Protocol differs from the implemented C2 norm, bounds, or scientific scope')


def hashes():
    paths = [HERE/n for n in ('README.md', 'protocol.json', 'requirements.txt', 'bounds.py',
                              'jets.py', 'review.py', 'audit.py', 'test_audit.py')]
    for directory, names in {
        'axis_core_attachment': ('README.md', 'bounds.py', 'evidence.json'),
        'axis_matching_audit': ('README.md', 'bounds.py', 'moments.py', 'annulus.py', 'evidence.json'),
        'axial_stress_audit': ('bounds.py',),
        'joined_stress_construction': ('README.md', 'construction.py', 'evidence.json'),
    }.items():
        paths += [RESULTS/directory/n for n in names]
    return {str(path.relative_to(RESULTS.parent)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def inherited_provenance():
    checks = {}
    for directory in ('axis_core_attachment', 'axis_matching_audit', 'joined_stress_construction'):
        data = json.loads((RESULTS/directory/'evidence.json').read_text())
        record = data.get('source_hashes', data.get('provenance'))
        if not isinstance(record, dict) or not record:
            raise ValueError('Missing inherited provenance: '+directory)
        checks[directory] = all(hashlib.sha256((RESULTS.parent/path).read_bytes()).hexdigest() == digest
                                for path, digest in record.items())
    return checks


def encode(v):
    if hasattr(v, '_mpi_'):
        return dict(lower=mp.nstr(lower(v), 120), upper=mp.nstr(upper(v), 120))
    if isinstance(v, mp.mpf):
        return mp.nstr(v, 90)
    if isinstance(v, Fraction):
        return str(v)
    if is_dataclass(v):
        return encode(asdict(v))
    if isinstance(v, dict):
        return {key: encode(value) for key, value in v.items()}
    if isinstance(v, (list, tuple)):
        return [encode(value) for value in v]
    return v


def generate():
    mp.mp.dps = 80
    protocol = json.loads((HERE/'protocol.json').read_text())
    validate_protocol(protocol)
    core = [core_transfer(d) for d in protocol['interval_digits']]
    rows = [normalization_constants(d) for d in protocol['interval_digits']]
    match = [matching_transfer(d) for d in protocol['interval_digits']]
    reviews = dict(normalization=normalization_review(), cancellation=cancellation_review(),
                   exponential=exponential_review(), axis_integrability=axis_integrability_review())
    provenance = inherited_provenance()
    gates = dict(actual_core_C2_and_pressure_integral_bounds=all(all(x['checks'].values()) for x in core),
                 five_normalized_C2_moment_bounds=all(all(x['checks'].values()) for x in rows),
                 exact_C2_matching_and_curvature=all(all(x['checks'].values()) for x in match),
                 all_derivative_controls=all(all(x['checks'].values()) for x in reviews.values()),
                 frozen_provenance=all(provenance.values()),
                 only_first_previous_obligation_closed=set(RESOLVED)|set(REMAINING) == set(previous.OPEN_OBLIGATIONS)
                    and not set(RESOLVED)&set(REMAINING),
                 actual_frequency_still_unverified=not FLAGS['actual_joined_profile_frequency_selected'])
    return dict(status=STATUS, protocol=protocol, flags=FLAGS, source_hashes=hashes(),
                core_transfer=core, normalized_rows=rows, matching=match, reviews=reviews,
                inherited_provenance=provenance, resolved_obligations=RESOLVED,
                remaining_obligations=REMAINING, gates=gates)


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS:
        raise ValueError('This checkpoint bounds only the actual C2 transfer, not the full stress or PDE')
    if record.get('resolved_obligations') != RESOLVED or record.get('remaining_obligations') != REMAINING:
        raise ValueError('The record silently changed the unresolved scientific scope')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A transfer or review check failed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-record', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = encode(generate())
    validate(result)
    if args.verify_record and json.loads(args.verify_record.read_text()) != result:
        raise SystemExit('Frozen evidence differs; no record written')
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(STATUS)
    print('Actual normalized target C2 <= 2.86e-17; matching coefficient curvature < 5.721e-14.')
    print(f"{len(result['gates'])} audit gates passed. Full stress, frequency, PDE and blow-up remain unverified.")


if __name__ == '__main__':
    main()
