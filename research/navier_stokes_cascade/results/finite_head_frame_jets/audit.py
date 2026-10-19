#!/usr/bin/env python3
"""Finite-head third-jet Borel and conditional radial spectral-frame C2 bounds."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from frame_jets import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='finite-head-third-jet-Borel-and-conditional-radial-spectral-frame-C2'
DEPENDENCIES=('collar_third_coefficient','cutoff_higher_jet_interface','regional_weighted_gate',
    'background_jet_foundation','current_frame_freezing','first_seminorm_bound')
FLAGS=dict(explicit_finite_head_geometric_tail_bound_from_original_diagonal_quantifiers=True,
    second_coefficient_third_norm_retained_unless_stronger_same_schedule_budget_verified=True,
    inverse_radius_full_shear_and_first_coefficient_product_jets_retained=True,
    conditional_current_native_spectral_C2_and_physical_radial_scale_bounds_derived=True,
    original_finite_head_proof_all_jet_Borel_requirements_pressure_and_construction_unchanged=True,
    historical_frozen_evidence_and_original_comparators_preserved=True,
    first_coefficient_plus_half_controls_every_derivative_order=False,
    actual_H23_or_enhanced_same_schedule_n2_third_budget_verified=False,
    abstract_scalar_profile_controls_are_actual_PDE_coefficients_or_residual_counterexamples=False,
    unchanged_cutoff_shape_jets_and_retained_scale_acceptance_verified=False,
    all_phase_normal_Cartesian_axial_time_or_generated_nonlinear_frame_gates_closed=False,
    complete_covariance_axial_pressure_means_and_residual_improvement_proved=False,
    independent_unconditional_core_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['bound_actual_second_coefficient_H23_or_verify_extra_same_schedule_n2_third_jet_budget_with_full_pressure_moment_sources',
    'supply_unchanged_cutoff_shape_jets_valid_total_current_jet_inputs_and_actual_value_mixed_retained_scale_acceptance',
    'join_all_phase_normal_potential_curl_Cartesian_axial_time_pulse_jets_geometry_supports_and_nonlinear_frame_inputs',
    'complete_full_wave_covariance_and_axial_velocity_pressure_mean_corrections',
    'prove_complete_viscous_cutoff_interaction_residual_and_independently_check_retained_foundations']
def protocol():
    return dict(parent_commit='79e54e2be7c0b250516e1e0cc2eb9462ffb3bf98',
        diagonal='original Cn,m q^(nh)<=2^-n for m<=n on whole support;chi independentX,eta;not a claim every derivative norm is budgeted at every n',
        finite_head='J=max(2,m),head n=1..J-1;error<=sum Hn,m z^(2n)+(z/2)^J/(1-z/2),z=q^h<=1;whole q2h constant<=sum head+2^(1-J)',
        third='B3(z)=H1,3+z^2 H2,3+z/[8(1-z/2)];H2,3 remains required unless SAME n2 third budget separately verified;old half tail still valid at orders1,2',
        scalar_control='auxiliary U2,N=a N^-2 psi(X)sin(N(X-1/2)),a=1/(1+2P1+P2);jets<=2 bounded1,central third=aN;profile-budget countercontrol,not actual PDE family',
        profile='E1=Rprof F1,deltaF=e/Rprof,g=(2X FX,Rprof UX);all radius product and lower shear terms retained;positive fixed native annulus',
        spectral='given current floors f<=F<=Fmax,-gtheta>=t,|g|<=G,D=|g|^2+2F gtheta>=d>0;N=g/|g|,K=rot N,c=-sqrt(D/(-2F gtheta));explicit quotient and normal C2 estimates from TOTAL jet inputs',
        physical='fixed physical z,t,X=r^2/(2q);radial factors q^-1/2 Rhi H1 and q^-1(H1+Rhi^2 H2);N,K,c band invariant,lambda additionally scales s^-1-h',
        pressure='actual H1,3,H2,3 depend on complete globally finalized coefficients and full lower-order pressure/moment memory;no velocity-only inference for relative stress or covariance',
        diagnostics='exact radius/shear,normal,quotient,physical-chain identities;five manufactured spectral profiles;four abstract second-coefficient controls;not selected enormous fields',
        scope='quantitative refinement of original finite-head proof,not an original Borel omission or schedule replacement;actual H23,cutoff/scale,all pulse maps,complete covariance/residual and unconditional foundations remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('finite_head_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','frame_jets.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py',
        'cutoff_higher_jet_interface/jet_bounds.py','collar_third_coefficient/third_bounds.py',
        'current_frame_freezing/mixed_bounds.py','first_seminorm_bound/seminorms.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Finite-head frame record changed')
    if not all(expected['gates'].values()):raise ValueError('Finite-head frame check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual H23, cutoff and complete pulse acceptance remain unverified.')
if __name__=='__main__':main()
