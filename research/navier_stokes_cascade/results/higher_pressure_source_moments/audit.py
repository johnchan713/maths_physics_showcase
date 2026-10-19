#!/usr/bin/env python3
"""Ordered full-support pressure moments and the conditional second repair."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from moments import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='ordered-full-support-pressure-source-moments-and-conditional-second-repair-bound'
DEPENDENCIES=('finite_head_frame_jets','signed_pressure_source_moments','coefficient_support_foundation',
    'global_axial_angular_bound','global_axial_source_envelope','background_summation')
FLAGS=dict(ordered_general_source_flux_and_signed_full_support_moments_derived=True,
    previous_order_axial_viscosity_starts_in_pressure_order_two_and_is_retained=True,
    finite_endpoint_radial_and_axial_boundary_terms_and_axis_trace_distinguished=True,
    explicit_conditional_Ceta1_moment_caps_use_global_angular_not_radial_jets=True,
    second_repair_retains_lower_F1_square_U1_square_and_U1F1_product=True,
    unchanged_linear_inverse_800_amplitude_slope_cost_and_full_pressure_memory_retained=True,
    historical_frozen_evidence_construction_cutoff_schedule_and_gauges_unchanged=True,
    actual_global_U1_Ceta3_and_U0_Ceta4_inputs_quantified=False,
    actual_second_inner_J3_H23_or_stronger_same_schedule_budget_verified=False,
    signed_moment_bounds_are_absolute_or_pointwise_source_bounds=False,
    weighted_radial_viscosity_cancellation_erases_shifted_axial_viscosity=False,
    actual_cutoff_jets_scale_acceptance_all_physical_pulse_and_support_bounds_closed=False,
    complete_covariance_axial_pressure_means_and_residual_improvement_proved=False,
    independent_unconditional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['quantify_actual_full_support_U1_Ceta3_U0_Ceta4_and_entire_axis_first_coefficient_inputs_for_second_inner_source',
    'derive_second_inner_J3_join_unchanged_cutoff_and_exact_repair_to_actual_H23_or_verify_same_schedule_extra_budget',
    'supply_actual_cutoff_jets_retained_value_scale_acceptance_and_all_physical_potential_curl_pulse_derivative_support_bounds',
    'complete_full_wave_covariance_axial_velocity_pressure_means_and_viscous_cutoff_interaction_residual',
    'independently_check_retained_core_regional_and_other_analytical_foundations']
def protocol():
    return dict(parent_commit='39dca194ea7dbe3f8e505ecaaa42729ca5e64754',
        hierarchy='k=n-1;lambda_j=2hj;all ordered i+j=k;w_k has -nu Z_(lambda_(k-1)-1-D) Z_(lambda_(k-1)-1) v_(k-1) for k>=1',
        flux='J_k=(V_k-2eta sum U_i V_j)/L+X sum v_i v_j-2nu V_kX;remainder=(-lambda_k v_k+Deta v_keta+d B_keta+2(lambda_k-A)eta B_k)/L+C_k/2-shifted viscosity',
        angular_operator='T_b=(d partial_eta+2b eta)/L;int X^j Z_(a-D)Z_a f=T_(a+j+1-D)T_(a+j+1) int X^j f only after full compact support',
        boundaries='axis 2nu v_k(0) retained;actual positive-order zero U axis jets imply it vanishes;finite R requires J_k(R),R J_k(R)+2nu V_k(R) and differentiated axial endpoint flux',
        support='same full Xv;functional M_j infinity zero;common endpoint independenteta;no Xcut I3 I4 truncation or incoming pressure reset',
        norms='Ceta1 signed totals from all U_i Ceta2,v_i Ceta2 and previous v Ceta3;regular average sufficient U_i Ceta3 and previous U Ceta4;not absolute or pointwise bounds',
        repair='Pi2_X=2F0F2+F1^2-w1/2;ordinary moments retain U1F1,U1^2,F1^2;unchanged exact B inverse800 with amplitude and lambda_patch inverse costs',
        diagnostics='48 Fraction compact C4 endpoint fixtures at k=0,1,3;192 integrated scalar jet identities;27 independent original floating hierarchy samples;not actual PDE profiles or interval certificates',
        scope='conditional global source/repair input only;actual U1 Ceta3,U0 Ceta4,second inner J3,H23,cutoff,scale,complete pulses,covariance,means,residual and independent foundations remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('higher_pressure_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','moments.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py',
        'signed_pressure_source_moments/moments.py','background_summation/hierarchy.py',
        'physical_residual_budget/residual.py','coefficient_support_foundation/support_bounds.py',
        'matching_root_foundation/matching_bounds.py','root_foundation/root_bounds.py')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
@lru_cache(None)
def reviews():return review()
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g=dict(r['checks']);g['historical_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Higher pressure moment record changed')
    if not all(expected['gates'].values()):raise ValueError('Higher pressure moment check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual global angular inputs and H23 remain unverified.')
if __name__=='__main__':main()
