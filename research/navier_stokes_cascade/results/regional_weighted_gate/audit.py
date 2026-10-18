#!/usr/bin/env python3
"""Finite-seminorm regional weighted gate, not actual coefficient norm evaluation."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from weighted_gate import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='finite-seminorm-weighted-input-and-regional-background-gap-join'
DEPENDENCIES=('background_gap_certificate','weighted_background_edges','background_jet_foundation',
 'weighted_stress_foundation','outer_edge_stability','coefficient_support_foundation','current_residual_ledger_audit')
FLAGS=dict(lumped_Cval_derived_from_actual_shared_budget_and_first_coefficient_seminorms=True,
 whole_small_q_higher_order_tail_bounded_by_half_q_power=True,
 common_velocity_and_relative_stress_cutoff_schedule_required=True,
 componentwise_to_Euclidean_norm_padding_and_all_four_Q_chart_factors_retained=True,
 exact_radial_inverse_radius_derivative_kept=True,
 compact_value_gate_joined_to_already_proved_nonnegative_outer_cone_monotonicity=True,
 large_outer_relative_error_not_mistaken_for_actual_cone_failure=True,
 full_radial_cutoff_curl_pressure_moments_Pi1_Pi2_and_regular_spatial_cutoff_terms_preserved=True,
 one_fixed_value_domain_all_fixed_jets_with_order_dependent_constants=True,
 earlier_heat_sign_and_canonical_seed_results_not_counted_as_new=True,
 historical_source_bytes_separate_from_fresh_diagnostic_replays=True,
 actual_first_coefficient_seminorms_and_huge_cutoff_scales_numerically_evaluated=False,
 all_primitive_core_regional_root_and_gap_estimates_independently_closed_here=False,
 velocity_only_cutoff_budget_implies_relative_stress_budget=False,
 outer_sign_theorem_rederived_from_scalar_controls=False,
 compact_constants_uniform_in_Xb_h_N_width_or_jet_order=False,
 full_domain_of_all_pulse_root_and_residual_maps_enlarged_here=False,
 actual_background_shear_frame_replaces_frozen_pulse_frame=False,
 historical_first_band_field_or_old_force_recertified=False,
 unconditional_pulse_complete_residual_single_field_force_or_blowup_proved=False,formal_certificate=False)
REMAINING=['derive_justified_actual_first_coefficient_and_compact_leading_gap_bounds_from_retained_core_regional_maps_and_check_the_combined_gate',
 'keep_historical_first_band_and_old_force_unverified_and_join_any_distinct_later_realization_with_its_own_controlled_schedules',
 'do_not_recount_completed_heat_sign_or_canonical_seed_results_or_promote_regional_value_bounds_to_unconditional_residual_or_blowup']

def protocol():
    return dict(parent_commit='0c7ff52961d301880d0078d462ea8fe402328855',
      common_schedule='one c_n dominates velocity and Euclidean relative-stress seminorms; all previous physical potential and higher-jet requirements retained',
      tail='sum_(n>=2) 2^-n q^(nh)=q^(2h)/[4(1-q^h/2)]<=q^(2h)/2 for 0<q<=1',
      first_data='A_E,A_EX,A_UX,A_W are actual first-coefficient sup seminorms on compact Kc; each B=A+1/2',
      pullback='dF=s^-A e/R; dgtheta=s^-A[(R/s)e_X-e/R]; dgz=s^-A(R/s)u_X; s^-1-h<5',
      constant='C_K=max(1,5B_E/Lr,5[Hr(B_EX+B_UX)+B_E/Lr],B_W)',
      regional_gate='compact C_Kc*q^(2h)<=r_Kc and exact outer nonnegative angular T1 preserve this cone interface; other global magnitude or map-domain gates remain',
      outer_limit='W1(Xb,+/-1)=2Xb/(1-2h)^2 can fail a sufficient global small-error test without actual cone failure',
      chart='whole q<=4Q requires ell>=2+log2(C_Kc/r_Kc)/(2h), combined with every prior accepted-tail inequality',
      memory='full pressure/moment sources, radial cutoff commutator, Pi1/Pi2,A0/A1,viscosity and all interactions unchanged',
      scope='conditional finite-seminorm input and regional cover only; actual seminorm values, primitive regional hypotheses and historical first band remain unverified')

@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('regional_weighted_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','weighted_gate.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py',HERE.parent/'background_jet_foundation/foundation.py',
        HERE.parent/'weighted_background_edges/edges.py',HERE.parent/'weighted_stress_foundation/edge_bounds.py',
        HERE.parent/'background_gap_certificate/gap_bounds.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
@lru_cache(None)
def reviews():return review()
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g=dict(r['checks'])
    g['historical_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Regional weighted record changed')
    if not all(expected['gates'].values()):raise ValueError('Regional weighted gate failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates; shared Borel budget, explicit chart constant and regional cone cover")
    print('Actual first-coefficient seminorm values, primitive regional hypotheses and unconditional theorem remain unpromoted.')
if __name__=='__main__':main()
