#!/usr/bin/env python3
"""Reproduce the loop estimate while preserving the unproved repair and PDE steps."""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp
from bounds import HERE, PROJECT, lower, upper, scalar_checks, derivative_ledger, error_ledger
from review import root_review, phase_review, modulation_review

STATUS = 'fixed-phase-loop-and-modulation-errors-bounded'
RESOLVED = ['actual_modulation_error_constants']
INHERITED = ['actual_target_C2_transfer', 'actual_compact_jet_envelope']
REMAINING = {'actual_repair_error_constants':
    'Bound the correction-to-state map by Ccorr on a positive-field neighborhood and justify its coefficient tolerance.'}
FLAGS = dict(actual_fixed_phase_loop_first_derivatives_bounded=True,
             actual_modulation_error_constants_bounded=True,
             relies_on_preceding_analytic_profile_and_envelope=True,
             actual_repair_error_constants_bounded=False,
             actual_joined_profile_frequency_selected=False,
             full_admissible_stress_realized=False,
             full_PDE_corrections_verified=False,
             smooth_force_verified=False, blowup_verified=False,
             formal_proof_assistant_certificate=False,
             independent_peer_review_completed=False)


def expected_protocol():
    return dict(status=STATUS, parent_commit='bb72b484feb0437800ea0a8cf9a2c7e1aefab527',
        source_pdf_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        input_envelope='A=C^4096', hierarchy=dict(G='exp(A^32)',H='exp(G^4)'),
        slow_derivative_order=1, moment_error_angular_order=1, pressure_state_error_order=0,
        Cstate='H^16', D_including_inverse_lambda='H^16',
        interval_digits=[80,110], diagnostic_digits=70,
        diagnostic_steps=['1e-5','1e-7','1e-9'], actual_joined_profile_frequency_selected=False)


def provenance():
    """Compare bytes, not the success flags of prior audits."""
    data = json.loads((HERE.parent/'compact_jet_envelope'/'evidence.json').read_text())
    old = data['source_hashes']
    if not old:
        raise ValueError('The inherited input record has no source hashes')
    return {p:hashlib.sha256((PROJECT/p).read_bytes()).hexdigest() == digest
            for p,digest in old.items()}


def hashes():
    names = ('README.md','protocol.json','requirements.txt','bounds.py','review.py','audit.py','test_audit.py')
    paths = [HERE/n for n in names]
    for directory in ('compact_jet_envelope','stress_realization_audit','joined_stress_construction'):
        paths += [HERE.parent/directory/n for n in ('README.md','evidence.json')]
    paths += [HERE.parent/'stress_realization_audit'/'loop.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def encode(v):
    if hasattr(v,'_mpi_'):
        return dict(lower=mp.nstr(lower(v),120),upper=mp.nstr(upper(v),120))
    if isinstance(v,mp.mpf):
        return mp.nstr(v,65)
    if isinstance(v,dict):
        return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):
        return [encode(x) for x in v]
    return v


def generate():
    protocol = json.loads((HERE/'protocol.json').read_text())
    if protocol != expected_protocol():
        raise ValueError('Derivative orders, scales, or scientific scope changed')
    scalars = [scalar_checks(d) for d in protocol['interval_digits']]
    reviews = dict(root=root_review(70),phase=phase_review(70),modulation=modulation_review(70))
    frozen = provenance()
    previous = json.loads((HERE.parent/'compact_jet_envelope'/'evidence.json').read_text())
    ledger = error_ledger()
    gates = dict(
        scalar_bounds_at_both_precisions=all(all(r['checks'].values()) for r in scalars),
        root_derivatives_and_zero_limits=all(reviews['root']['checks'].values()),
        moving_inverse_phase_and_omission_control=all(reviews['phase']['checks'].values()),
        exact_modulation_and_failure_controls=all(reviews['modulation']['checks'].values()),
        state_and_moment_error_ledger=all(ledger['checks'].values()),
        inherited_source_bytes_unchanged=all(frozen.values()),
        previous_scope_preserved=set(previous['remaining_obligations']) == set(RESOLVED)|set(REMAINING),
        all_diagnostics_labelled_manufactured=all(r['manufactured'] for r in reviews.values()))
    return dict(status=STATUS,flags=FLAGS,protocol=protocol,source_hashes=hashes(),
                inherited_resolved_obligations=INHERITED,resolved_obligations=RESOLVED,
                remaining_obligations=REMAINING,scalar_bounds=scalars,
                derivative_ledger=derivative_ledger(),error_ledger=ledger,
                reviews=reviews,inherited_provenance=frozen,gates=gates)


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS:
        raise ValueError('Loop error bounds do not certify the repair, frequency, or PDE')
    if (record.get('resolved_obligations') != RESOLVED
            or record.get('inherited_resolved_obligations') != INHERITED
            or record.get('remaining_obligations') != REMAINING):
        raise ValueError('An unresolved obligation was removed or silently promoted')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A bound, diagnostic, or provenance check failed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-record',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = encode(generate())
    validate(result)
    if args.verify_record and json.loads(args.verify_record.read_text()) != result:
        raise SystemExit('Evidence differs; no output written')
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(STATUS)
    print('Cstate and transformed D bounded; inverse-phase derivative and 1/lambda retained.')
    print('Repair neighborhood, accepted frequency, full PDE, smooth force and blow-up remain unverified.')


if __name__ == '__main__':
    main()
