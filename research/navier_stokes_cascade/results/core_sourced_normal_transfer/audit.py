#!/usr/bin/env python3
"""Freeze original normal-floor transfer and full-pressure denominator bounds."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_normal_length_closure','core_phase_coefficient_budget','core_first_sourced_inverse_budget')
STATUS='original-normal-floor-transferred-to-sourced-inverse-without-stronger-phase-beta-gate'
FLAGS={
    'original_inverse_normal_square_4_over_bmin_squared_reaches_sourced_inverse':True,
    'independent_stronger_phase_beta_gate_removed_for_this_original_tail_route':True,
    'all_order_positive_denominator_degrees_velocity_4j_plus_2_pressure_potential_4j_plus_3':True,
    'uniform_denominator_change_adds_constants_not_scale_powers_at_fixed_order':True,
    'unprojected_pressure_source_moving_normal_and_cutoff_products_retained':True,
    'new_analytical_assumptions_added':False,
    'actual_input_enclosures_or_unknown_historical_band_certified':False,
    'tangent_propagator_complete_covariance_mean_or_residual_closed':False,
    'historical_fields_bands_lengths_supports_or_evidence_changed':False,
    'unconditional_pulses_or_blowup_proved':False,
}
REMAINING=[
    'retain valid original primitive bounds and existing original accepted-tail membership; no stronger phase-beta gate is needed for this sourced inverse route',
    'establish homogeneous tangent propagator Gaussian-ratio bound on the SAME canonical entrance full path and source domains',
    'retain actual full background phase and source jet bounds and full parameter-path support; fixed-input denominator comparison alone does not bound feedback source dependence',
    'join generated complete covariance fresh axial velocity and full-memory pressure means, then all physical scales viscosity cutoff and interaction residual estimates',
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
    paths += [HERE.parent/n/'budgets.py' for n in DEPENDENCIES]
    paths += [HERE.parent/'frame_column_foundation/frame_budget.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_original_floor_degree_pressure_and_transfer_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   three_historical_records_bound=len(bound)==3,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Original normal-floor transfer denominator degree pressure or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} inherited source checks")


if __name__=='__main__':main()
