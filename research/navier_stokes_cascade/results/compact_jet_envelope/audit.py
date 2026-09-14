#!/usr/bin/env python3
"""Reproduce the compact input envelope without promoting the proposed frequency."""
import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import mpmath as mp
from bounds import (HERE, RESULTS, c2, previous, scalar_bounds, margin_bounds,
                    exponent_ledger, lower, upper)
from review import source_and_jet_review, cutoff_review, shear_sign_review

STATUS = 'actual-compact-input-envelope-bounded'
FLAGS = dict(actual_target_C2_transfer_bounded=True,
             actual_matching_coefficients_C2_bounded=True,
             actual_compact_jet_envelope_verified=True,
             relies_on_preceding_analytic_core_argument=True,
             actual_modulation_error_constants_bounded=False,
             actual_repair_error_constants_bounded=False,
             actual_joined_profile_frequency_selected=False,
             full_admissible_stress_realized=False,
             full_PDE_corrections_verified=False,
             smooth_force_verified=False, blowup_verified=False,
             formal_proof_assistant_certificate=False,
             independent_peer_review_completed=False)
INHERITED_RESOLVED = ['actual_target_C2_transfer']
RESOLVED = ['actual_compact_jet_envelope']
REMAINING = {k:v for k,v in previous.OPEN_OBLIGATIONS.items()
             if k not in INHERITED_RESOLVED+RESOLVED}
DEPENDENCIES = ('axis_core_attachment','axis_matching_audit','angular_c2_transfer',
                'joined_stress_construction','axial_stress_audit','intermediate_decay_audit')


def validate_protocol(protocol):
    expected = dict(status=STATUS, parent_commit='028333027590ae2bad71f5be74d699f4cb979c3d',
        source_pdf_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_pages=[26,27,28,29,30,31,145,146,147,148,150,154,155,156],
        scale='Md64', interval_digits=[80,110],
        radial_derivative='d/dy = X*d/dX, y=log(X)',
        mixed_norm='sum_{r+s<=k} sup(abs(d_y^r d_eta^s f))/(r!s!)',
        field_and_moment_order=2, stress_coordinate_order=1, envelope='C^4096',
        actual_joined_profile_frequency_selected=False)
    if protocol != expected:
        raise ValueError('Protocol changed the derivative orders, domain, envelope or scientific scope')


def inherited_provenance():
    checks = {}
    for directory in DEPENDENCIES:
        data = json.loads((RESULTS/directory/'evidence.json').read_text())
        record = data.get('source_hashes',data.get('provenance'))
        if not isinstance(record,dict) or not record:
            raise ValueError('Missing frozen source hashes: '+directory)
        def source(path):
            # Older audits used paths relative to results, newer ones to the project.
            return (RESULTS.parent if path.startswith('results/') else RESULTS)/path
        checks[directory] = all(hashlib.sha256(source(path).read_bytes()).hexdigest() == digest
                                for path,digest in record.items())
    return checks


def hashes():
    paths = [HERE/n for n in ('README.md','protocol.json','requirements.txt','bounds.py',
                             'jets.py','review.py','audit.py','test_audit.py')]
    for directory in DEPENDENCIES:
        paths += [RESULTS/directory/n for n in ('README.md','evidence.json')]
        paths += [v for v in (RESULTS/directory/'bounds.py',RESULTS/directory/'construction.py')
                  if v.exists()]
    paths += [RESULTS/'stress_realization_audit'/'README.md']
    return {str(path.relative_to(RESULTS.parent)):hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def encode(v):
    if hasattr(v,'_mpi_'):
        return dict(lower=mp.nstr(lower(v),120),upper=mp.nstr(upper(v),120))
    if isinstance(v,mp.mpf):
        return mp.nstr(v,90)
    if isinstance(v,Fraction):
        return str(v)
    if is_dataclass(v):
        return encode(asdict(v))
    if isinstance(v,dict):
        return {key:encode(value) for key,value in v.items()}
    if isinstance(v,(list,tuple)):
        return [encode(value) for value in v]
    return v


def generate():
    mp.mp.dps = 80
    protocol = json.loads((HERE/'protocol.json').read_text())
    validate_protocol(protocol)
    scalars = [scalar_bounds(d) for d in protocol['interval_digits']]
    margins = [margin_bounds(d) for d in protocol['interval_digits']]
    ledger = exponent_ledger()
    inherited_C2 = [dict(core=c2.core_transfer(d)['checks'],
                        matching=c2.matching_transfer(d)['checks'])
                    for d in protocol['interval_digits']]
    cutoff_constants = dict(second=c2.matching.certificate(110)['checks']['second_step_derivative_below_20000'],
                            third=previous.derivative_constants(110)['checks']['step_third_below_1e6'])
    reviews = dict(source_and_jets=source_and_jet_review(), cutoff=cutoff_review(), shear_sign=shear_sign_review())
    provenance = inherited_provenance()
    gates = dict(actual_regional_scalar_bounds=all(all(r['checks'].values()) for r in scalars),
                 actual_original_margin_bounds=all(all(r['checks'].values()) for r in margins),
                 field_moment_and_pressure_exponents=all(ledger['checks'].values()),
                 preceding_actual_C2_bounds=all(all(v.values()) for r in inherited_C2 for v in r.values()),
                 cutoff_derivative_constants=all(cutoff_constants.values()),
                 independent_equation_paths_and_failure_controls=all(all(r['checks'].values()) for r in reviews.values()),
                 frozen_provenance=all(provenance.values()),
                 only_compact_input_newly_closed=set(INHERITED_RESOLVED+RESOLVED)|set(REMAINING)
                    == set(previous.OPEN_OBLIGATIONS) and len(REMAINING) == 2,
                 actual_frequency_and_full_construction_unverified=not any(FLAGS[k] for k in
                    ('actual_joined_profile_frequency_selected','full_admissible_stress_realized',
                     'full_PDE_corrections_verified','smooth_force_verified','blowup_verified')))
    return dict(status=STATUS,protocol=protocol,flags=FLAGS,source_hashes=hashes(),
                scalar_bounds=scalars,margin_bounds=margins,exponent_ledger=ledger,
                preceding_C2_checks=inherited_C2,cutoff_constants=cutoff_constants,
                reviews=reviews,inherited_provenance=provenance,
                inherited_resolved_obligations=INHERITED_RESOLVED,resolved_obligations=RESOLVED,
                remaining_obligations=REMAINING,gates=gates)


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS:
        raise ValueError('The compact input envelope does not certify the loop, frequency or PDE')
    if (record.get('resolved_obligations') != RESOLVED
        or record.get('inherited_resolved_obligations') != INHERITED_RESOLVED
        or record.get('remaining_obligations') != REMAINING):
        raise ValueError('The unresolved modulation or repair obligation was changed')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('An analytic constant, equation review or provenance gate failed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-record',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = encode(generate())
    validate(result)
    if args.verify_record and json.loads(args.verify_record.read_text()) != result:
        raise SystemExit('Frozen evidence differs; no record written')
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(STATUS)
    print('Actual compact input envelope A=C^4096; narrow-cutoff losses retained.')
    print(f"{len(result['gates'])} audit gates passed. Modulation, repair errors and frequency remain open.")


if __name__ == '__main__':
    main()
