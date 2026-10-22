#!/usr/bin/env python3
"""Freeze the original-frame physical tangent propagator and inverse join."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_sourced_normal_transfer','core_normal_length_closure','frame_column_foundation','sourced_inverse_foundation')
STATUS='original-frame-budgets-close-physical-Gaussian-ratio-tangent-propagator-and-length-inputs'
FLAGS={
    'original_CE_plus_Cd_supplies_integrated_logarithmic_norm_bound':True,
    'physical_Euclidean_endpoint_frame_conversion_proved':True,
    'uniform_all_nonzero_harmonic_value_bound_keeps_extra_actual_damping':True,
    'tangent_constant_and_original_L_le_S_inserted_in_unchanged_inverse':True,
    'canonical_frame_mismatch_not_an_independent_zero_entry_inverse_obstruction':True,
    'full_pressure_moving_normal_and_original_denominator_retained':True,
    'new_analytical_assumptions_or_acceptance_gates_added':False,
    'actual_primitive_enclosures_or_unknown_historical_band_certified':False,
    'homogeneous_primary_frames_or_covariance_roots_identified':False,
    'complete_covariance_fresh_means_source_support_or_residual_closed':False,
    'historical_fields_bands_lengths_supports_or_evidence_changed':False,
    'unconditional_pulses_or_blowup_proved':False,
}
REMAINING=[
    'retain actual validity and uniformity of ORIGINAL full background/frame primitive bounds and same accepted-tail/domain; do not reopen T or L upper bound as independent premises on this route',
    'justify actual coefficient/source jet and whole-parameter-path support bounds including generated complete covariance',
    'carry fresh axial velocity full Pi1/Pi2 pressure memory and moving moments; resolve homogeneous primary-frame/root identification where those identities require it',
    'complete all physical-scale derivative support viscosity cutoff and interaction residual estimates before unconditional claims',
]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    paths=[HERE/n for n in ('README.md','budgets.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_sourced_normal_transfer/budgets.py',HERE.parent/'core_normal_length_closure/budgets.py',
              HERE.parent/'frame_column_foundation/frame_budget.py',HERE.parent/'frame_column_foundation/frame_identity.py',
              HERE.parent/'core_first_sourced_inverse_budget/budgets.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_energy_frame_gauge_harmonic_and_inverse_join_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   four_historical_records_bound=len(bound)==4,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Physical tangent propagator frame energy inverse join or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} inherited source checks")


if __name__=='__main__':main()
