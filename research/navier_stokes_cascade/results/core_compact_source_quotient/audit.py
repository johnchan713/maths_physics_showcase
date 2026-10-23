#!/usr/bin/env python3
"""Freeze exact compact-source identities, pressure memory and scoped A1 gain."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_compact_input_sharpness','core_next_sourced_mean_tail_join','compact_mean_foundation',
              'compact_mean_corrections','core_fast_inverse_budget','core_coupled_defect_budget')
STATUS='exact-compact-source-quotient-and-conditional-actual-axial-potential-A1-gain'
FLAGS={
    'original_full_drift_boundary_and_compact_projection_identities_derived':True,
    'current_radial_source_reduced_but_K0_pressure_memory_restored':True,
    'same_actual_compact_axial_potential_retains_gd_minus_A1_gd':True,
    'time_and_both_viscosities_annihilated_only_in_A1_image':True,
    'whole_axial_nonlinear_covariance_pressure_flux_and_moving_cutoff_kept':True,
    'conditional_A1_error_B_plus_one_two_kappa_above_previous_source_route':True,
    'direct_normalized_A1_estimate_uses_no_nonstationary_integration':True,
    'postwave_A1_and_final_A0_stages_and_original_6q_12q_domains_kept_distinct':True,
    'radial_swirl_obstruction_respects_divergence_and_both_Haar_velocity_moments':True,
    'radial_source_actual_high_jet_support_and_pressure_inputs_validated':False,
    'whole_Ez_pressure_or_complete_regular_E_order_promoted':False,
    'next_generic_signed_cycle_or_controlled_infinite_diagonal_completed':False,
    'canonical_source_or_actual_residual_lower_bound_produced':False,
    'absolute_decay_unconditional_pulses_or_blowup_proved':False,
    'historical_evidence_actual_fields_h_width_roots_carriers_or_band_changed':False
}
REMAINING=[
    'validate actual r0/H1 coarse and high-jet weighted/profile norms, exact full potential/stress extension and support identities on every original band with all current covariance and full pressure memory',
    'resolve the surviving actual centripetal A0 and compensated pressure source bounds; divergence and the two velocity moments alone do not remove them',
    'advance the next generic signed cycle on the enlarged domains, retaining older A0/root/normal request DAGs and the full additive tail ledger',
    'obtain sufficient complete regular orders and controlled infinite diagonal, global localization/forcing regularity and velocity growth before any blow-up claim'
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
    paths += [HERE.parent/'core_compact_input_sharpness/bounds.py',HERE.parent/'core_next_sourced_mean_tail_join/budgets.py',
              HERE.parent/'core_signed_prefix_source_budget/budgets.py',HERE.parent/'core_fast_inverse_budget/proof_checks.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {k:canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_drift_projection_pressure_and_conditional_A1_controls_pass=bool(r['checks']) and all(r['checks'].values()),
            six_historical_records_bound=len(bound)==6,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Compact quotient, full pressure memory, actual potential route, scopes or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
