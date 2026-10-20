#!/usr/bin/env python3
"""Native complete pressure J4 and a field-only fourth slow-jet boundary."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from boundary import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='native-pressure-joint-four-and-explicit-fourth-slow-frequency-boundary'
DEPENDENCIES=('higher_jet_budget','actual_fifth_native','mixed_shape_cone_bridge',
    'fourth_loop_shape_gate','fourth_pressure_input_boundary','fifth_moment_volterra',
    'premodulation_mixed_shear','fourth_target_normalization','repair_state_bounds',
    'cutoff_higher_jet_interface','modulation_repair_all_order')
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
FLAGS=dict(complete_PREMOD_native_pressure_jointJ4_A1547_derived=True,
    original_five_densities_full_axis_Cp_axial_pressure_retained=True,
    exact_normalized_fourth_frequency_and_all_exponential_partitions_retained=True,
    missing_slow_fourth_Gamma_enters_field_majorant_linearly=True,
    eight_remaining_fifth_mixed_field_inputs_identified=True,
    fixed_N_functional_nonimplication_not_actual_profile_counterexample=True,
    same_root_widths_supports_frequency_and_historical_evidence_unchanged=True,
    actual_cutoff_Z_or_Gamma_inequality_verified=False,
    actual_fifth_mixed_shear_inputs_enclosed=False,
    unconditional_FINAL_compact_field_or_global_leading_C4_verified=False,
    first_C3_second_H23_uniform_actual_pulses_complete_covariance_means_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['bound_actual_fifth_mixed_field_indices_2_3_3_2_4_1_5_0_for_E_and_U_then_full_slow_fourth_loop_Gamma',
    'enclose_unchanged_actual_cutoff_Z_and_Gamma_field_gate;full_FINAL_pressure_flux_jets_still_require_their_own_inputs',
    'prove_FINAL_global_leading_C4_first_C3_second_H23_uniform_all_scale_actual_pulses_supports_complete_covariance_axial_pressure_means_and_full_viscous_cutoff_interaction_residual']
def protocol():
    return dict(parent_commit='dca2a2cdbadd9f3e0a56f3a49fc0d4abf6012441',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        native='same PREMOD rectangle,A=C4096>=2^16,Bn=A256;field_fullPi_inverseJ4 and W_p1_p2J3<=Bn;original value floors and native/core hypotheses retained',
        pressure='differentiate original forward map using all five density flows;full Pi=Pi_axis+Cp;radial-containing normalized pressure rows<A771,pure angularC4<=A1546,sum<A1547<G;NOT FINAL pressure',
        shear='fourth native shears a=1-2Ey/E,bs=2Uy/E have principal coefficients -2(r+1)/E0,+2(r+1)/E0;remaining E,U fifth mixed indices (2,3),(3,2),(4,1),(5,0) not enclosed',
        frequency='same N=1+floorH32;(P_N)_(r,s)=sum_c N^(c-1)P_(r-c,s,c);under sufficient Z parent full P,g=Pphi jointJ3<=H and angular fourth<=H;Gamma is sum of pure-slow fourth rows r>=1 of both primitives',
        field='same leading five-row repair Ceta4<=2000H16/N;unit-mass width.1 bumps charge radialJ4<1e106<Bn;fieldJ4<=32Bn N3 H+2Bn Gamma/N+2;K=H64 unchanged,separate fieldK2=H128 sufficient if Gamma<=N H128/(4Bn);actual Z,Gamma unverified',
        control='fixed protected chi and ANY fixed finite N;zero-mean B_n=n^-3 chi cos(n(y-y0))sin(2pi phi);P,g jointJ3<=1e211<H,pure angular fourth zero,normalized composed fourth/n tends1/(24N);functional nonimplication NOT actual C9 or global pressure-compatible counterexample',
        scope='conditional-native analytical pressure estimate,conditional field-only transfer and precise unverified slow-jet boundary;finite symbolic,rational,interval and manufactured checks are not continuum proof;no FINAL global derivative,pulse,covariance,mean,residual or blowup acceptance')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(None)
def history():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('native_mixed_fourth_historical_bytes',p/'audit.py')
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
    files=[HERE/n for n in ('README.md','boundary.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Native mixed fourth boundary record changed')
    if not all(expected['gates'].values()):raise ValueError('Native mixed fourth boundary check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual cutoff and Gamma gates and FINAL/pulse acceptance remain open.')
if __name__=='__main__':main()
