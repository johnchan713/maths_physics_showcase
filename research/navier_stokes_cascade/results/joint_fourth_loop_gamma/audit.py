#!/usr/bin/env python3
"""Full joint fourth loop and Gamma under an explicit unverified target premise."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from joint import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='joint-fourth-moving-loop-Gamma-and-conditional-compact-field-gate'
DEPENDENCIES=('native_fifth_mixed_cover','native_mixed_fourth_boundary','fourth_loop_shape_gate',
    'mixed_shape_cone_bridge','fourth_target_normalization','modulation_repair_all_order',
    'repair_state_bounds','regional_cone_join','global_moment_stress',
    'actual_cutoff_pullback','cutoff_composition_obstruction','compact_pressure_bridge')
HISTORICAL_RECORDS=(
    'activated_collar_foundation',
    'actual_band_pulse_bounds',
    'actual_cutoff_pullback',
    'actual_fifth_native',
    'angular_c2_transfer',
    'annular_exterior_foundation',
    'annular_first_correction',
    'auxiliary_geometry_audit',
    'axial_root_foundation',
    'axial_stress_audit',
    'axis_core_attachment',
    'axis_first_coefficient',
    'axis_matching_audit',
    'axis_second_coefficient',
    'background_gap_certificate',
    'background_jet_foundation',
    'background_summation',
    'canonical_comparison_join',
    'canonical_endpoint_join',
    'canonical_finite_cycle_join',
    'canonical_force_join',
    'canonical_initial_mean_join',
    'canonical_local_realization',
    'canonical_primary_foundation',
    'canonical_value_domain',
    'certified_tail_reconstruction',
    'coefficient_support_audit',
    'coefficient_support_foundation',
    'collar_operator_envelope',
    'collar_third_coefficient',
    'compact_jet_envelope',
    'compact_mean_corrections',
    'compact_mean_foundation',
    'compact_pressure_bridge',
    'comparison_foundation',
    'complete_covariance_initial_means',
    'complex_mean_relaxation',
    'complex_signed_covariance_update',
    'core_all_order_foundation',
    'core_contraction_audit',
    'core_growth_activation',
    'core_pressure_identification',
    'current_complex_foundation',
    'current_endpoint_jets',
    'current_frame_freezing',
    'current_localized_force',
    'current_potential_realization',
    'current_residual_ledger_audit',
    'current_solution_comparison',
    'current_target_cycle_bridge',
    'current_target_moment_audit',
    'cutoff_composition_obstruction',
    'cutoff_cone_scale_gate',
    'cutoff_higher_jet_interface',
    'cutoff_repair_derivative_budget',
    'cutoff_shape_feasibility',
    'endpoint_foundation',
    'fifth_moment_volterra',
    'finite_head_frame_jets',
    'finite_induction_foundation',
    'finite_stage_induction',
    'first_correction_construction',
    'first_seminorm_bound',
    'first_sourced_foundation',
    'first_sourced_wave_update',
    'force_foundation',
    'fourth_angular_axial_interface',
    'fourth_loop_shape_gate',
    'fourth_pressure_input_boundary',
    'fourth_target_normalization',
    'frame_column_foundation',
    'full_potential_realization',
    'gaussian_mass_foundation',
    'geometric_pulse_replay',
    'global_axial_angular_bound',
    'global_axial_source_envelope',
    'global_moment_stress',
    'global_pressure_moment_foundation',
    'global_solution_comparison',
    'global_swirl_envelope',
    'harmonic_closure_foundation',
    'heat_exterior_audit',
    'higher_jet_budget',
    'higher_pressure_source_moments',
    'inner_collar_activation',
    'intermediate_decay_audit',
    'joined_stress_construction',
    'leading_cone_margin',
    'leading_frame_foundation',
    'leading_frame_freezing',
    'localized_force_extension',
    'loop_modulation_bounds',
    'matching_root_foundation',
    'mixed_shape_cone_bridge',
    'modulation_repair_all_order',
    'native_fifth_mixed_cover',
    'native_mixed_fourth_boundary',
    'nonlinear_axis_pilot',
    'outer_edge_stability',
    'outer_pressure_pilot',
    'paper_profile_audit',
    'physical_pulse_foundation',
    'physical_pulse_transfer',
    'physical_residual_budget',
    'physical_wave_interface',
    'post_pulse_stress_audit',
    'post_signed_mean_relaxation',
    'premodulation_mixed_shear',
    'prescribed_continuation_foundation',
    'primary_covariance_foundation',
    'pulse_moment_audit',
    'pulse_stress_audit',
    'real_signed_covariance_update',
    'realization_foundation',
    'reduced_continuation_margins',
    'reference_fifth_transfer',
    'regional_cone_join',
    'regional_stress_floor',
    'regional_weighted_gate',
    'regular_source_foundation',
    'repair_state_bounds',
    'residual_cycle_closure',
    'root_foundation',
    'same_tilt_deficit',
    'second_coefficient_source',
    'shape_join_pressure',
    'sharp_current_scale_floor',
    'signed_pressure_source_moments',
    'source_envelope_edge_audit',
    'sourced_inverse_foundation',
    'stress_realization_audit',
    'weighted_background_edges',
    'weighted_stress_foundation',
)
FLAGS=dict(complete_joint_fourth_signed_root_moving_source_and_primitive_estimates_derived=True,
    full_three_external_variable_inverse_and_primitive_operator_three_retained=True,
    Gamma_less_than_2H_CONDITIONAL_on_joint_target_J4_bound_derived=True,
    same_frequency_exact_repair_compact_FIELD_J4_less_than_K_squared_CONDITIONAL=True,
    original_five_pressure_memories_axis_axial_and_U_DeltaW_retained=True,
    old_K_N_roots_widths_supports_cutoffs_and_historical_evidence_unchanged=True,
    actual_Z_or_joint_target_J4_condition_verified=False,
    unconditional_actual_Gamma_or_FINAL_compact_field_accepted=False,
    complete_FINAL_fifth_angular_moment_or_global_leading_flux_pressure_C4_verified=False,
    first_C3_second_H23_uniform_actual_pulses_complete_covariance_means_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['enclose_unchanged_actual_cutoff_Z_or_directly_prove_joint_target_J4_bound;all_new_Gamma_and_field_conclusions_remain_conditional',
    'supply_FINAL_fifth_angular_moment_inputs_and_full_leading_global_flux_pressure_C4;PREMOD_J5_is_not_FINAL_J5',
    'prove_first_C3_second_H23_all_scale_actual_pulses_supports_complete_covariance_axial_pressure_means_and_full_viscous_cutoff_interaction_residual']
def protocol():
    return dict(parent_commit='4bca4deafaf306be4fbd87e4d33a99ef7bfa5707',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        native='same PREMOD full pressure p2J4<A1547,native EJ5_UJ5<A5,a_bsJ4<A262,inverse_aJ4<A1054,tsJ4<A1316,vsJ4<A1579<G;not old Bn pressure label or FINAL jets',
        condition='J:s joint(y,eta)J4<=G;parent raw shape bound1e12 Z deltaL^-8 A13685;old sufficient explicit Z gate power20673 implies J with spare6988logA;actual Z and J unverified',
        scalar='retain q5,R4 and regular signed F(mu,p2),Fpartials4<exp64G,Fmu inverse<exp10G;no division by p2,mu or Vmu at zero;two slow colors counts4,20,96,448,mu exponents80,240,400,560',
        moving='full y_eta_theta fixed angle sourceJ4<exp20000G,wJ4<exp10000G,PhiJ4<exp10010G;three external colors y_eta_varphi counts6,45,324,2268 and inverse exponents11000,33000,55000,77000;positive ThetaJ4<exp80000G,moving sourceJ4<exp350000G',
        primitive='exact leading zero means with full variance and axial density;phase-free integration norm2,positive phase normalized derivative g_(k-1)/k norm1;FULL joint operator<=3,primitiveJ4<exp400000G<H;Gamma<2H CONDITIONAL on J',
        field='same N=1+floorH32 and exact five-row root;fieldJ4<=32Bn N3 H+2Bn Gamma/N+2<H128 under J;Gamma is not dropped,old K=H64 unchanged,K2 is separate fourth FIELD budget;repair derivative width.1 and all pressure memories retained',
        boundary='original radial flux normalized fourth angular principal term -5(1-eta2)/(1-2h eta2)*M_(0,5);FINAL fifth angular memory and full global flux_pressureC4 remain unbounded;no PREMOD-to-FINAL promotion',
        scope='conditional analytical joint fourth Gamma and field gate;finite symbolic/rational/interval diagnostics not continuum proof;no actual cutoff or unconditional FINAL global pulse covariance mean residual or blowup acceptance')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(None)
def history():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('joint_fourth_gamma_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
        h=m.verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
    finally:sys.path[:]=prior
    bindings=h['bindings']
    if not bindings or not all(bindings.values()):raise ValueError('Historical source binding failed')
    if not HISTORICAL_RECORDS:raise ValueError('Fixed parent manifest required')
    records={str(p.relative_to(PROJECT)):sha(p) for p in (HERE.parent/n/'evidence.json' for n in HISTORICAL_RECORDS)}
    return dict(all_bindings_verified=True,binding_count=len(bindings),
        binding_result_sha256=hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest(),frozen_records=records)
def sources():
    files=[HERE/n for n in ('README.md','joint.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py',
        'native_fifth_mixed_cover/mixed_bounds.py','native_mixed_fourth_boundary/boundary.py',
        'fourth_loop_shape_gate/loop_fourth.py','mixed_shape_cone_bridge/bridge.py')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in files}
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(review());h=history();g=dict(r['checks']);g['historical_source_bytes']=h['all_bindings_verified']
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Joint fourth loop Gamma record changed')
    if not all(expected['gates'].values()):raise ValueError('Joint fourth loop Gamma check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual cutoff and unconditional FINAL/global/pulse acceptance remain open.')
if __name__=='__main__':main()
