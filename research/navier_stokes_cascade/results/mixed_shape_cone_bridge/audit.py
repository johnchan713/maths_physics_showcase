#!/usr/bin/env python3
"""Audit fresh joint loop and same-root cone bounds under an explicit shape gate."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from bridge import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='mixed-third-loop-and-exact-repair-cone-bridge-under-explicit-shape-gate'
DEPENDENCIES=('actual_cutoff_pullback','fourth_loop_shape_gate','higher_jet_budget','premodulation_mixed_shear',
    'compact_jet_envelope','loop_modulation_bounds','repair_state_bounds','regional_cone_join',
    'compact_pressure_bridge','global_moment_stress')
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
    'modulation_repair_all_order',
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
FLAGS=dict(mixed_target_bound_includes_actual_shape_factor_and_width_cost=True,
    sufficient_shape4_gate_implies_new_mixed3_gate=True,
    joint_third_moving_loop_rederived_without_borrowing_unbudgeted_joint_H=True,
    complete_phase_primitive_operator_and_three_variable_inverse_retained=True,
    same_N_exact_five_row_root_and_conditional_compact_value_cone_join_derived=True,
    full_axis_Cp_axial_and_U_DeltaW_pressure_memory_retained=True,
    roots_supports_width_frequency_and_historical_evidence_unchanged=True,
    actual_cutoff_budget_or_M_or_Z_gate_verified=False,
    unconditional_actual_compact_cone_or_global_leading_admissibility_verified=False,
    FINAL_full_mixed_leading_C4_first_C3_second_H23_or_uniform_all_scale_pulses_verified=False,
    complete_covariance_axial_pressure_means_residual_gain_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['enclose_unchanged_actual_Z_or_complete_actual_mixed3_and_angular4_target_composition;neither_gate_currently_verified',
    'prove_FINAL_full_mixed_leading_C4_first_C3_second_H23_and_retained_scale_acceptance_without_unbudgeted_cutoff_inputs',
    'prove_uniform_actual_physical_pulses_common_supports_complete_wave_covariance_axial_pressure_means_and_full_viscous_cutoff_interaction_residual']
def protocol():
    return dict(parent_commit='66ae6450fb32aad2720462f45ca29c47eb7cd80f',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        native='same compact PREMOD rectangle,native jointJ3 inputs<=Bn=A256 from nativeJ4 field/moment ledger;sharp lower-order/value A unchanged;full p2 pressure memory,not p1 as a root argument',
        mixed_shape='actual sJ3<=1e6 Z3 A2308 deltaL^-6;M:logZ3+6log10+2308logA+6A20<=A32 is sufficient,not verified',
        stronger_gate='prior Z4 gate logZ4+12log10+20673logA+8A20<=A32 implies M since Z3<=Z4 and gap=6log10+18365logA+2A20>0;pure angular S4 or necessary pullback gate alone not promoted',
        moving='fresh factorial colored-jet collapse;mu degrees80,240,400;Theta degrees11000,33000,55000 in y_eta_varphi;full primitiveJ3<exp400000G<H with full phase operator<=3',
        state='fixed-phase derivatives then Ny evaluation retain fast radial term and exp(-Acal/N);Cstate=Dtarget=H16;Ccorr=A128;ballA^-128;coefficient<=2000H16/N;total state error<G^-64/2 under M',
        cone='same N=1+floorH32;all compact gaps>=G^-4/2 and a>=G^-3/2 CONDITIONAL on M and retained regional/core/inverse premises;protected collars and exact five-row memory join to compensated exterior',
        memory='same actual axis pressure,Cp and axial velocity;all five exact matching functions retain derivatives;U_DeltaW outside zero-U repair patch is not dropped;no positive-order future-pressure reset',
        scope='conditional analytical value/cone bridge with original maps and exact manufactured diagnostics;no actual shape acceptance,FINAL full mixed leadingC4,firstC3,secondH23,unconditional physical pulse covariance residual or blowup')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(None)
def history():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('mixed_shape_bridge_historical_bytes',p/'audit.py')
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
    files=[HERE/n for n in ('README.md','bridge.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','higher_jet_budget/higher_bounds.py',
        'fourth_loop_shape_gate/loop_fourth.py','repair_state_bounds/bounds.py')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in files}
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(review());h=history();g=dict(r['checks']);g['historical_source_bytes']=h['all_bindings_verified']
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Mixed shape cone bridge record changed')
    if not all(expected['gates'].values()):raise ValueError('Mixed shape cone bridge check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual shape gate and unconditional compact/FINAL/pulse acceptance remain open.')
if __name__=='__main__':main()
