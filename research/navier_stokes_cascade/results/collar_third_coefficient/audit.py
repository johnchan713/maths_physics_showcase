#!/usr/bin/env python3
"""Activated third leading jets and explicit inner J3; no actual cutoff acceptance."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from third_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='actual-activated-third-leading-jet-and-conditional-inner-coefficient-J3'
DEPENDENCIES=('cutoff_higher_jet_interface','collar_operator_envelope',
    'activated_collar_foundation','core_all_order_foundation','first_seminorm_bound')
FLAGS=dict(actual_activated_leading_X3_eta2_envelope_derived_relative_to_retained_core=True,
    full_matrix_and_pressure_source_xi_derivative_envelopes_derived=True,
    explicit_inner_J3_on_positive_cutoff_collar_from_same_sparse_Picard_solution=True,
    angular_Cauchy_losses_and_singular_diagonal_radial_derivatives_charged=True,
    full_pressure_shifted_viscosity_future_moments_and_same_repair_interface_retained=True,
    historical_frozen_evidence_and_construction_choices_unchanged=True,
    radial_analyticity_across_activation_or_axis_inverse_radius_estimate_used=False,
    unchanged_cutoff_L3_and_actual_retained_scale_acceptance_verified=False,
    illustrative_C_M_inputs_are_selected_actual_core_values=False,
    finite_difference_diagnostics_are_interval_or_continuum_certificates=False,
    all_higher_coefficient_orders_or_complete_physical_pulse_jets_closed=False,
    complete_covariance_axial_pressure_means_and_residual_improvement_proved=False,
    independent_unconditional_core_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['supply_unchanged_cutoff_shape_jets_and_full_global_repair_inverse_inputs_then_check_actual_mixed_retained_scale_gate',
    'propagate_every_required_higher_order_and_verify_all_physical_scale_Cartesian_axial_time_pulse_jets_geometry_and_supports',
    'retain_all_jet_Borel_potential_curl_cutoff_costs_and_nonlinear_frame_acceptance_inputs',
    'complete_full_wave_covariance_and_axial_velocity_pressure_mean_corrections',
    'prove_complete_viscous_cutoff_interaction_residual_and_independently_check_retained_foundations']
def protocol():
    return dict(parent_commit='2fcf54abdc2b3faf661d73e0feae79197f96c565',
        leading='same actual squared-flat activation t=C^-20,tc=t/(10sqrt(log C)),ell=tc/4<t/20;core radial orders<=3,zero-free reciprocal,actual Fa=Fr exp J,Ua_y=(1-ea)Ur_y',
        third_activation='Jyyy<=30400S/t^2+128000S^2/t+50000S^3;full Bell and logarithmic Jacobian terms;no radial analyticity across flat join',
        separated='actual K3,H3 bounded by rational factors times E<=2^ceil(50ell S);positive operator value/xi derivative polynomials degree<=2 in H3,so both factor bounds multiply by E^2',
        operator='positive value/xi derivative pairs of every original A0,A1 row and source;regular w0 keeps radial viscosity,full shifted axial operators keep coefficient derivatives;nu=1',
        pressure='source=(0,0,0,-xi*w0,-2CF,2eta X*w0/L-2CU);both pressure columns retained;local estimate never erases full future moment/inverse costs through Xv',
        Picard='C>=1,Cxi>=0 on2delta;z=2Cb,beta=e z^2/(delta/2),M=(z+beta)exp beta on3delta/2,using retained cross-point sparse derivative block proof',
        derivative='B=C+3/a,Bxi=Cxi+3/a^2,M1=2M/delta,M2=8M/delta^2;V1=BM+CM1+C;V1eta=(C/delta)M+(B+C/delta)M1+CM2+C/delta;V2=Bxi M+B V1+Cxi M1+C V1eta+Cxi',
        J3='Ceta=1+2/delta;J3=Ceta[V2/(8a^3)+3V1/(8a^4)+3M/(8a^5)],0<a<=sqrt(Xa) on cutoff collar;all bounds affine in unevaluated M',
        diagnostic='18 original-system manufactured matrix/source value and centered-difference xi derivative comparisons,exact product/coordinate checks,toy M<=2^4000 at C=2,b=2,loss1/8;not selected solution',
        scope='inner J3 input discharged relative to retained core/Picard;unchanged cutoff L3,actual scale,all pulses,higher coefficients,complete covariance/residual and unconditional foundations remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('third_collar_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','third_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py',
        'first_correction_construction/system.py','collar_operator_envelope/envelope.py',
        'cutoff_higher_jet_interface/jet_bounds.py','physical_residual_budget/review.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Third collar record changed')
    if not all(expected['gates'].values()):raise ValueError('Third collar check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual cutoff/scale and complete pulses remain unverified.')
if __name__=='__main__':main()
