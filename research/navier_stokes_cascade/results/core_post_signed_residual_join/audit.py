#!/usr/bin/env python3
"""Freeze the conditional current full-source and physical-request tail join."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_post_signed_mean_closure','core_generic_signed_covariance_join','core_multiharmonic_mean_residual_budget',
              'core_fast_row_gain_join','core_physical_tail_join','regular_source_foundation',
              'current_residual_ledger_audit','core_signed_target_edge_budget')
STATUS='conditional-post-signed-current-nonzero-source-and-full-physical-request-tail-join'
FLAGS={
    'new_current_nonzero_source_derived_from_actual_signed_and_fresh_mean_differences':True,
    'old_regular_source_added_not_cancelled_by_homogeneous_signed_update':True,
    'full_signed_pressure_retains_moving_normal_and_K_with_original_projected_solve':True,
    'complete_uncut_curl_viscosity_spatial_cutoffs_and_both_transport_orientations_kept':True,
    'actual_A1_velocities_distinct_from_current_additive_compact_residuals':True,
    'regular_nonzero_D_and_full_regular_mean_C0_give_conditional_quantitative_seed':True,
    'physical_compact_request_graph_not_replaced_by_normalized_low_output_graph':True,
    'old_DAGs_raw_kernel_floors_and_all_high_derivative_requests_remain_deferred':True,
    'separate_oriented_primitive_constant_restored_in_common_sum_shift':True,
    'archived_lower_shifts_retained_without_assuming_pole_monotonicity':True,
    'same_band_physical_flatness_recipes_and_assigned_base_tail_retained':True,
    'new_h_same_band_gate_physical_field_or_width_selected':False,
    'actual_weighted_high_jet_support_background_or_assigned_base_bounds_validated':False,
    'nonzero_source_promoted_to_mean_order_C0':False,
    'physical_absolute_decay_or_residual_lower_bound_proved':False,
    'full_induction_infinite_field_unconditional_pulses_or_blowup_proved':False,
    'historical_frozen_evidence_changed':False,
}
REMAINING=[
    'validate actual same-family background/root/coefficient/profile/support/endpoint and assigned base/exterior inequalities through all new and older nested finite requests',
    'apply the NEXT constrained sourced inverse to this exact current supported regular source and derive the full next cycle on common expanded domains',
    'prove higher complete regular orders satisfying h*alpha>N+Lambda+2m+eta; current D still does not certify physical absolute decay',
    'join a controlled infinite diagonal, localization product jets, smooth forcing and full velocity growth before any blow-up claim',
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
    paths += [HERE.parent/'core_post_signed_mean_closure/budgets.py',HERE.parent/'core_generic_signed_covariance_join/budgets.py',
              HERE.parent/'core_multiharmonic_mean_residual_budget/budgets.py',HERE.parent/'core_fast_row_gain_join/budgets.py',
              HERE.parent/'core_physical_tail_join/bounds.py',HERE.parent/'core_signed_target_edge_budget/budgets.py',
              HERE.parent/'finite_induction_foundation/finite_checks.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}

def canonical(value):
    if isinstance(value,Fraction):return str(value)
    if isinstance(value,dict):return {k:canonical(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [canonical(v) for v in value]
    return value

def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_current_pressure_Fourier_source_and_physical_request_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   eight_historical_records_bound=len(bound)==8,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Current full source, signed pressure, physical domains, limits or provenance changed')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    payload=json.dumps(r,indent=2,sort_keys=True)+'\n'
    with a.output.open('x') as f:f.write(payload)
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")

if __name__=='__main__':main()
