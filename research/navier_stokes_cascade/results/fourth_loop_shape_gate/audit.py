#!/usr/bin/env python3
"""Audit fourth scalar/moving loop and explicitly conditional cutoff-shape gate."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from loop_fourth import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='fourth-angular-signed-root-moving-phase-and-conditional-shape-gate'
DEPENDENCIES=('actual_fifth_native','premodulation_mixed_shear','fourth_pressure_input_boundary',
    'fourth_target_normalization','fourth_angular_axial_interface','higher_jet_budget',
    'loop_modulation_bounds','modulation_repair_all_order','repair_state_bounds')
FLAGS=dict(joint_q5_removable_R4_and_signed_F_partial4_budgets_derived=True,
    fourth_signed_root_and_full_moving_phase_inverse_angular_budgets_derived_conditionally_on_actual_s_C4_le_G=True,
    angular_C4_primitives_below_unchanged_H_under_shape_condition_derived_without_old_joint_C3_budget=True,
    same_full_five_row_target_normalizers_inverse_lambda_and_exact_small_root_retained=True,
    sufficient_conditional_angular_compact_C4_gate_has_headroom_at_same_N_and_K=True,
    full_native_pressure_axis_Cp_and_axial_memory_enters_p2_without_reset=True,
    precise_unverified_variance_shape_inequality_kept_separate_from_radial_sigma=True,
    roots_supports_width_frequency_pressure_gauges_and_historical_evidence_unchanged=True,
    actual_variance_cutoff_shape_enclosed_or_s_C4_condition_verified=False,
    unconditional_actual_compact_angular_gate_or_radial_state_cone_acceptance_verified=False,
    FINAL_full_mixed_leading_C4_first_C3_second_H23_or_uniform_all_scale_pulses_closed=False,
    full_covariance_axial_pressure_means_viscous_cutoff_interaction_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['identify_actual_unchanged_variance_cutoff_Z_or_independent_s_C4_le_G_composition_proof;plateaus_smoothness_and_radial_sigma_do_not_supply_it',
    'join_remaining_radial_state_cone_and_FINAL_full_leading_first_coefficient_and_second_H23_inputs_and_retained_scale_gates',
    'prove_all_physical_pulse_derivatives_common_supports_complete_wave_covariance_axial_pressure_means_full_viscous_cutoff_interaction_residual_and_independent_foundations']
def protocol():
    return dict(parent_commit='66ff8c49955a1b17e04a617d42681633cc19b51b',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        domain='same compact PREMOD native data and actual circle loop,eta[-1,1],uniform y as parameter;new estimates angular and auxiliary-phase only,not fourth y or full joint J4',
        parameters='A=C4096>=2^16,G=expA32,H=expG4,N=1+floorH32,d0=1/(16A),deltaL=exp(-A20);all unchanged;nativep2C4<=A1546 and nativea_bs_ts_vsC4<=A2583<G',
        unproved_condition='S:actual sCeta4<=G;Z sufficient:logZ+12log10+20673logA+8A20<=A32;Z=max1,raw derivatives1..4 of unchanged variance shape on unit transition;Z1e100 auxiliary only',
        scalar='q joint5<exp4G,R joint4<exp5G,R2C4<exp10G and >=exp(-2G)/2,sqrtR2C4<exp60G,signedFpartial4<exp64G,Fmu^-1<exp10G;no division by p mu s or Vmu',
        implicit='normalized mu exponents80,240,400,560;phase exponents11000,33000,55000,77000;all partitions and p2 fourth jet retained;not p1 as a separate root argument',
        moving='fixed-angle t joint4<exp4000G,w<exp10000G,Phi<exp10010G,Theta positive C4<exp80000G;full moving g sources<exp350000G;periodic primitives angularC4<exp400000G<H CONDITIONALLY on S',
        compact='all five original target normalizers<=A5,raw moment affine coefficients<=A515,ZH Za Zb<=A522,alpha beta ch<=A525<G;C4 target<=3A522H/N<H16/N gives SAME small root and r;Kknown<A259<G;conditional compactcost<G+1<H<K=H64',
        memory='actual native p2 includes all five entering fifth moments and same axis datum;only modulation differences cancel identical axis;exact repair restores Cp as well as other rows;no positive-order future-pressure reset',
        scope='conditional analytic estimate plus exact finite and manufactured numerical diagnostics;actual shape S/Z remains unverified;no unconditional compact/radial/FINAL/global/pulse/covariance/residual/blowup or independent foundation claim')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('fourth_loop_shape_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','loop_fourth.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','higher_jet_budget/higher_jets.py',
        'fourth_target_normalization/target_bounds.py','fourth_angular_axial_interface/angular_bounds.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Fourth loop shape record changed')
    if not all(expected['gates'].values()):raise ValueError('Fourth loop shape check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual variance shape/S condition and unconditional compact/FINAL acceptance remain open.')
if __name__=='__main__':main()
