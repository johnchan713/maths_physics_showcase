#!/usr/bin/env python3
"""Audit an actual native crossing and necessary, not sufficient, shape estimate."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from pullback import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='actual-axial-taper-crossing-and-necessary-cutoff-pullback-budget'
DEPENDENCIES=('cutoff_composition_obstruction','fourth_loop_shape_gate','premodulation_mixed_shear',
    'annular_exterior_foundation','actual_fifth_native','higher_jet_budget')
HISTORICAL_RECORDS=(
    'activated_collar_foundation',
    'actual_band_pulse_bounds',
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
FLAGS=dict(actual_native_axial_taper_crosses_entire_variance_window_non_degenerately=True,
    necessary_shape_budget_extracted_from_actual_target_with_complete_inverse_C4=True,
    actual_pressure_not_replaced_by_zero_pressure_diagnostic=True,
    full_pressure_axis_Cp_axial_memory_roots_width_supports_frequency_and_history_unchanged=True,
    retained_earlier_annular_match_and_Md64_parameter_hypotheses_explicit=True,
    actual_selected_Z_enclosed_or_S_condition_verified=False,
    necessary_shape_inequality_is_sufficient_for_S=False,
    actual_S_condition_or_primitive_H_bound_disproved=False,
    FINAL_global_jets_unconditional_pulses_full_covariance_means_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['enclose_unchanged_actual_Z_or_prove_complete_actual_S_le_G;actual_crossing_rules_out_angular_flatness_but_does_not_supply_Z',
    'join_radial_state_cone_FINAL_leading_C4_first_C3_second_H23_with_same_pressure_memory_roots_width_and_frequency',
    'prove_uniform_physical_pulse_derivatives_supports_complete_covariance_axial_pressure_means_and_full_viscous_cutoff_interaction_residual']
def protocol():
    return dict(parent_commit='8c34b682cd1d1ee593f190ae46efabe65a647f1f',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        native='actual PREMOD axial taper E=P1*f*exp(-y/2),U=k(y)*eta,Md64;earlier annular edits ended with exact native match;slice y0=exp32-1 before final eleven flat units',
        geometry='a=2,vs=2+beta*eta^2*(1+eta^2)^2,beta=.5exp(y0-64)/P1^2;exp(-65)A^-2<beta<1;A=C4096,deltaL=exp(-A20);r=sqrt(deltaL/beta)<1e-6',
        pullback='eta(t)=inverse_of_eta_plus_eta_cubed(r*sqrt((1+t)/8));actual s(eta(t))=sqrt(deltaL)*K(t)*sqrt(3-t)/4;full inverse through C4 retained',
        bound='normalized inverse coefficients<=1,3,16,105;etaC4<=126r<1;K_C4<=8deltaL^-1/2*S;Z<=192deltaL^-1/2*S',
        necessary='S<=G implies logZ<=A32+A20/2+log192;necessary only,not actual Z enclosure or acceptance of S or failure of H',
        sufficient_retained='prior logZ+12log10+20673logA+8A20<=A32 implies S<=G;not evaluated for actual cutoff',
        memory='actual full p2 never set to zero;no axis/Cp/axial datum,future pressure tail,moment root,support,width C^-20,cutoff or N altered',
        scope='analytical actual-native estimate relative to retained matching and parameter premises;80/110 outward comparisons plus manufactured noninterval inverse diagnostics;no unconditional FINAL physical pulse covariance residual or blowup')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(None)
def history():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('actual_cutoff_pullback_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
        h=m.verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
    finally:sys.path[:]=prior
    bindings=h['bindings']
    if not bindings or not all(bindings.values()):raise ValueError('Historical source binding failed')
    if not HISTORICAL_RECORDS:raise ValueError('Fixed historical manifest required')
    records={str(p.relative_to(PROJECT)):sha(p) for p in (HERE.parent/n/'evidence.json' for n in HISTORICAL_RECORDS)}
    return dict(all_bindings_verified=True,binding_count=len(bindings),
        binding_result_sha256=hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest(),frozen_records=records)
def sources():
    files=[HERE/n for n in ('README.md','pullback.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','higher_jet_budget/higher_bounds.py','axis_core_attachment/bounds.py')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in files}
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(review());h=history();g=dict(r['checks']);g['historical_source_bytes']=h['all_bindings_verified']
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Actual cutoff pullback record changed')
    if not all(expected['gates'].values()):raise ValueError('Actual cutoff pullback check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual Z/S and unconditional compact/FINAL acceptance remain open.')
if __name__=='__main__':main()
