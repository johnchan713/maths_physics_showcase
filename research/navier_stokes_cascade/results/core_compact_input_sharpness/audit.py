#!/usr/bin/env python3
"""Freeze the exact scalar compact-input obstruction and its limited scope."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_next_sourced_mean_tail_join','core_fast_inverse_budget','compact_mean_foundation','core_physical_tail_join')
STATUS='sharp-isotropic-finite-jet-compact-input-obstruction-not-actual-source-validation'
FLAGS={
    'original_compact_profiles_moving_cutoff_hull_and_drift_retained':True,
    'real_smooth_comparison_sources_have_zero_Haar_and_all_Haar_moments':True,
    'Pell_comparison_handles_any_large_original_rounded_M_without_retuning':True,
    'source_radial_phase_cancels_original_compact_shift_exactly':True,
    'uniform_scalar_M_minus_p_rule_requires_J_at_least_m_plus_two_p':True,
    'original_sufficient_raw_rule_is_three_above_power_boundary':True,
    'original_epsilon_flooring_gives_finite_jet_threshold_up_to_polynomials':True,
    'parent_normalized_and_physical_route_requests_are_three_to_five_above_threshold':True,
    'comparison_sources_are_actual_canonical_sources_or_uniform_all_jet_class_members':False,
    'physical_derivative_loss_two_m_proved_necessary_for_all_arguments':False,
    'equality_boundary_or_three_derivative_upper_gap_resolved':False,
    'current_source_high_jet_support_pressure_input_inequalities_validated':False,
    'subsequent_generic_signed_cycle_or_infinite_diagonal_completed':False,
    'actual_residual_lower_bound_absolute_decay_unconditional_pulses_or_blowup_proved':False,
    'historical_evidence_fields_profiles_width_h_carriers_or_accepted_band_changed':False
}
REMAINING=[
    'validate actual current Qcirc gr_final and Qcirc Epostwave_z high-jet weighted-to-raw norms on every original band and whole support hull, retaining all covariance, actual A1 velocities and full pressure memory',
    'if reducing those derivative requests, prove additional quantitative structure of the SAME canonical sources which excludes the scalar resonant comparison family',
    'advance the next generic signed correction and fresh means on the enlarged original operator domains; no later order is supplied here',
    'obtain sufficient complete regular orders and controlled infinite diagonal, localization/forcing regularity and full velocity growth before any blow-up claim'
]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks

def sources():
    paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_next_sourced_mean_tail_join/budgets.py',HERE.parent/'core_fast_inverse_budget/proof_checks.py',
              HERE.parent/'core_physical_tail_join/bounds.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}

def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {k:canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v

def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_operator_sharpness_and_parent_request_controls_pass=bool(r['checks']) and all(r['checks'].values()),
            four_historical_records_bound=len(bound)==4,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Scalar compact obstruction, exact requests, scope or provenance changed')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")

if __name__=='__main__':main()
