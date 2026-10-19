#!/usr/bin/env python3
"""Global axial/source inputs; no full swirl or residual acceptance claim."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from envelopes import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='global-axial-pointwise-source-and-shifted-viscosity-envelope'
DEPENDENCIES=('cutoff_repair_derivative_budget','global_axial_angular_bound','higher_jet_budget',
    'collar_operator_envelope','signed_pressure_source_moments','first_seminorm_bound')
FLAGS=dict(same_compact_loop_J3_supplies_mixed_third_field_bound_without_fourth_loop_input=True,
    selective_axial_log_jets_bounded_with_unchanged_positive_N_loss=True,
    complete_reference_axial_taper_pulse_and_both_affine_end_bumps_bounded=True,
    global_radial_axial_derivatives_joined_without_dividing_at_axis=True,
    exact_divergence_lowers_pointwise_source_demand_to_total_mixed_order_three=True,
    separate_pointwise_and_absolute_source_bounds_include_viscosity_and_angular_derivatives=True,
    complete_shifted_axial_U_operator_includes_differentiated_coordinate_coefficients=True,
    full_future_Xv_source_and_nonzero_reserved_plateau_pressure_memory_retained=True,
    all_roots_supports_parameters_width_frequency_gauges_and_both_mean_memories_unchanged=True,
    historical_frozen_evidence_preserved=True,
    actual_giant_parameter_profile_numerically_materialized=False,
    global_swirl_and_shifted_swirl_inputs_closed=False,
    actual_cutoff_shape_combined_Cstar_tilt_and_scale_gate_verified=False,
    retained_core_and_regional_premises_independently_proved_unconditional=False,
    unconditional_actual_pulses_complete_residual_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['bound_actual_global_swirl_F_radial_and_shifted_axial_inputs_through_later_reference_and_heat_regions',
    'identify_or_bound_unchanged_inner_coefficient_cutoff_shape_and_join_full_Cstar_to_same_tilt_and_scale_gates',
    'discharge_all_scale_actual_pulse_full_covariance_mean_and_viscous_interaction_residual_obligations',
    'independently_verify_retained_exact_core_and_regional_foundations_before_unconditional_claims']
def protocol():
    return dict(parent_commit='dfcf62462a6a10da6313be53638f3ab2b25d9583',
        compact='weighted J3_N pullback<=H/N; same leading repair J3<H^20/N; native E,U J3<H^4*N^2<H^70',
        selective='|U_y|2<H^2; |U_yy-U_y|1<H^3*N, after axis region; original N=1+floor(H^32)',
        reference='same taper and complete pulse U J3<1e110<H, using actual width .3 unnormalized end bumps and both affine coefficients',
        global_inputs='KU=max(Kaxis,2Ahi,352); R1=max(Kaxis,H^2/Xa),R2=max(Kaxis,H^3*N/Xa^2)',
        source='VX and VXX divergence identities lower derivative demand; separate w0 C1 and absolute integrals retain viscosity=1 and full future support Xv',
        axial='full shifted Z_(k-D)Z_k U bound 64Xv^2 R2+704Xv R1+1664KU; all coordinate derivatives retained',
        scope='conditional leading axial/source inputs; actual global swirl, cutoff shape, complete first-stress acceptance, full residual and blow-up remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('global_axial_source_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','envelopes.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','signed_pressure_source_moments/moments.py',
        'higher_jet_budget/higher_jets.py','pulse_moment_audit/bounds.py','outer_pressure_pilot/schedule.py','physical_residual_budget/residual.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Global axial/source record changed')
    if not all(expected['gates'].values()):raise ValueError('Global axial/source check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; global swirl, cutoff shape and full same-field stress acceptance remain unverified.')
if __name__=='__main__':main()
