#!/usr/bin/env python3
"""Audit a local moving-loop obstruction without promoting it to the actual profile."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from composition import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='local-loop-composition-obstruction-not-an-actual-profile-counterexample'
DEPENDENCIES=('fourth_loop_shape_gate','premodulation_mixed_shear','actual_fifth_native')
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
FLAGS=dict(local_native_C4_envelopes_do_not_imply_target_or_moving_primitive_C1_budget=True,
    full_moving_inverse_and_phase_mean_subtraction_retained_in_obstruction=True,
    exact_complete_fourth_composition_principal_and_cross_terms_checked=True,
    actual_pressure_memory_roots_supports_width_frequency_and_frozen_history_unchanged=True,
    parent_conditional_fourth_angular_gate_retained_without_scope_promotion=True,
    actual_retained_cutoff_or_composition_S_enclosed=False,
    local_test_family_is_actual_global_pressure_compatible_profile=False,
    actual_primitive_H_bound_disproved=False,
    unconditional_compact_FINAL_pulses_covariance_means_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['enclose_unchanged_actual_cutoff_Z_or_prove_complete_actual_s_C4_le_G_with_crossing_or_cancellation_control',
    'join_radial_state_cone_FINAL_leading_C4_first_C3_second_H23_at_same_roots_frequency_and_width',
    'prove_uniform_physical_pulses_supports_complete_covariance_axial_pressure_means_and_full_viscous_cutoff_interaction_residual']
def protocol():
    return dict(parent_commit='117160fffb9824fb26831b72ab8ce03638eba082',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        family='auxiliary local eta[-1,1],0<delta<=1,n>=1 integer;a=E=1,vs=2+3delta/16+delta*eta/8,bs=sqrt(vs-1),p2=0 scalar slice;not global pressure realization',
        cutoff='auxiliary zeta_n(v)=1-sigma(n*(8*(v-2)/delta-3/2)+1/2);permitted plateaus;never substituted for actual variance cutoff',
        bounds='native normalized C4 norms<3 independent of n;sfirst>=2sqrt5*n*sqrt(delta);moving zero-mean A_eta norm>sqrt10*n*sqrt(delta)/(6pi)',
        moving='signed p2=0 branch,t=ts+sqrt2*s*sin(theta),v=vs+s^2;fixed-origin inverse and zero-mean primitive;moving derivative=(v/vs)*frozen derivative at theta=pi',
        fourth='raw complete Deta4(zeta(vs)*h);principal h*zeta4*(vs_eta)^4;no claimed actual transversality or cancellation exclusion',
        memory='no actual pressure datum or five-moment/future tail changed;no physical field,root,cutoff,frequency,width or support edited',
        scope='analytical nonimplication for generic interface only;manufactured floating quadrature diagnostics;actual S/Z and actual primitive H neither proved nor disproved;parent implication still conditional')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(None)
def history():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('composition_obstruction_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
        h=m.verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
    finally:sys.path[:]=prior
    bindings=h['bindings']
    if not bindings or not all(bindings.values()):raise ValueError('Historical source binding failed')
    records={str(p.relative_to(PROJECT)):sha(p) for p in (HERE.parent/n/'evidence.json' for n in HISTORICAL_RECORDS)}
    return dict(all_bindings_verified=True,binding_count=len(bindings),
        binding_result_sha256=hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest(),frozen_records=records)
def sources():
    files=[HERE/n for n in ('README.md','composition.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in files}
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(review());h=history();g=dict(r['checks']);g['historical_source_bytes']=h['all_bindings_verified']
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Composition obstruction record changed')
    if not all(expected['gates'].values()):raise ValueError('Composition obstruction check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual S/Z and unconditional compact/FINAL acceptance remain open.')
if __name__=='__main__':main()
