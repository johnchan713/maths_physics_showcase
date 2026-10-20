#!/usr/bin/env python3
"""Tame FINAL fifth moments and explicit complete compact fourth-map input gate."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from fifth import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='tame-FINAL-fifth-moment-and-complete-compact-fourth-map-interface'
DEPENDENCIES=('joint_fourth_loop_gamma','native_fifth_mixed_cover','actual_fifth_native',
    'fifth_moment_volterra','fourth_target_normalization','modulation_repair_all_order',
    'repair_state_bounds','fourth_pressure_input_boundary','global_moment_stress',
    'compact_jet_envelope','native_mixed_fourth_boundary','mixed_shape_cone_bridge')
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
    'joint_fourth_loop_gamma',
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
FLAGS=dict(tame_FINAL_fifth_angular_field_moment_and_same_root_bounds_affine_in_Lambda5_derived=True,
    all_five_density_normalizer_and_incoming_axis_pressure_memories_retained=True,
    fifth_constant_Q_repair_uses_existing_pointwise_Jacobian_not_new_C5_contraction=True,
    complete_compact_joint_fourth_pressure_flux_and_ordinary_conversion_interface_derived=True,
    bounded_fourth_functional_inputs_do_not_automatically_control_partial_fifth_moments=True,
    original_roots_K_N_widths_supports_cutoffs_gauge_and_history_unchanged=True,
    actual_fifth_primitive_Lambda5_or_sufficient_F5_gate_enclosed=False,
    actual_Z_or_joint_target_J4_condition_verified=False,
    unconditional_COMPLETE_compact_or_global_leading_C4_accepted=False,
    first_C3_second_H23_uniform_actual_pulses_complete_covariance_means_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['enclose_actual_pulled_primitive_fifth_angular_Lambda5_including_all_new_fifth_native_cutoff_moving_inputs;F5_not_verified',
    'verify_unchanged_actual_joint_target_or_raw_cutoff_budget_and_join_compact_full_maps_to_remaining_global_regions_and_retained_scales',
    'prove_first_C3_second_H23_coefficient_reserved_mean_supports_all_scale_actual_pulses_complete_covariance_axial_pressure_means_and_full_viscous_cutoff_interaction_residual']
def protocol():
    return dict(parent_commit='827e2c1f9915d28a0246830b9992cdddb2d1a382',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        premises='retained native/core/regional/floor/root premises and parent J:joint targetJ4<=G;J and actual cutoff unverified;full primitiveJ4<H,FINAL fieldJ4<H128,rootCeta4<=2000H16/N,nativefieldJ5<A5,nativemoment_axisCeta5<A8',
        high='Lambda5=sum of sup normalized fifth eta jets of BOTH actual pulled primitives at fixed external phase;not bounded by H or assigned old J4 budget;eta evaluation phase Ny does not remove it',
        modulation='exact fifth exponential seven partitions,top high jet linear;Dmod=16Bn(H+Lambda5)/N;full fifth product difference<=6Bn D using lowerC4 smallness;all five moments including U2_E2_Cp retained',
        normalization='full incoming momentCeta5 error<A515(H+Lambda5)/N;reciprocal polynomial costs5,25 through5;normalize all5 rows before inverse_lambda M-J transform,operator<=4A6,targetCeta5<A522(H+Lambda5)/N',
        repair='SAME LATER leading constant B_Q root,pointwise Jacobian inverse<=2000;fifth recurrence4 lower pairs,known tail3.2e15 H32/N2;rootCeta5<(H20+A524Lambda5)/N;earlier Q_eta not frozen and no C5 contraction or N reselection',
        moments='FINAL differencefieldsCeta5<(H24+A784Lambda5)/N,lowerC4<H20/N;all moment_fullPi Ceta5<A8+(H26+A1044Lambda5)/N;radial rows use actual densityJ4<H258,full momentJ5<mathcalM=H260+A1044Lambda5/N;no axis reset',
        maps='FINAL inverseEJ4<H514,Hx inverseJ4<H516;complete W,Qs,Ns,p1,p2,V maps linear in mathcalM;common logJ4<H660mathcalM,full ordinary conversion<H4,ordinaryC4<H664mathcalM;F5 Lambda5<=N H260/A1044 sufficient for separate H1024 compact budget,not old K replacement or global acceptance',
        control='fixed finite N and positive-lobe protected chi;B_n=n^-4 chi(y)sin(n eta)sin2pi_phi has uniform joint4 but partial M5=c_N n/(120N),complete normalized fourth flux/n tends-c_N/(24N);functional not actual C9/cone/profile failure',
        scope='conditional analytical affine fifth-moment and compact complete fourth-map interface;actual Lambda5,F5,J,Z not verified;finite checks not continuum proof;no unconditional global physical pulse covariance mean residual or blowup')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(None)
def history():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('final_fifth_moment_historical_bytes',p/'audit.py')
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
    files=[HERE/n for n in ('README.md','fifth.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py',
        'joint_fourth_loop_gamma/joint.py','fifth_moment_volterra/volterra.py',
        'fourth_target_normalization/target_bounds.py','modulation_repair_all_order/README.md')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in files}
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(review());h=history();g=dict(r['checks']);g['historical_source_bytes']=h['all_bindings_verified']
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('FINAL fifth moment interface record changed')
    if not all(expected['gates'].values()):raise ValueError('FINAL fifth moment interface check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual Lambda5/cutoff and unconditional global/pulse acceptance remain open.')
if __name__=='__main__':main()
