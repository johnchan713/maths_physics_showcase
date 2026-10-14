#!/usr/bin/env python3
"""Conditional actual-band pulse/primary/linear bounds; full cycle stays open."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
import numpy as np
from diagnostics import (band_frame_review, cone_and_edge_review, derivative_cutoff_review,
                         fixture_first_band, linear_remainder_review, partition_scaling_review,
                         phase_identity_review, pulse_covariance_review, support_review,
                         threshold_and_barrier_review)

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'conditional-actual-band-pulse-primary-covariance-and-linear-bounds'
FLAGS = dict(conditional_actual_background_to_band_phase_frame_estimate_supplied=True,
             conditional_actual_pulse_positive_lower_bound_and_all_fixed_derivatives_supplied=True,
             conditional_primary_covariance_positive_inverse_and_edge_bounds_supplied=True,
             conditional_signed_inverse_at_fixed_primary_supplied=True,
             conditional_full_linear_wave_residual_gain_supplied=True,
             conditional_real_primary_curl_covariance_remainder_gain_supplied=True,
             source_harmonic_temporal_zero_extension_required_for_global_linear_gain=True,
             relies_on_inherited_continuum_background_leading_gaps_and_support_geometry=True,
             actual_global_background_numerically_evaluated=False,
             actual_profile_constants_or_actual_qstar_numerically_selected=False,
             global_auxiliary_color_centers_and_rectangles_materialized=False,
             independent_foundation_and_support_geometry_review_completed=False,
             actual_mean_correctors_constructed=False,
             full_nonlinear_interaction_and_residual_cycle_verified=False,
             full_physical_stress_wave_construction_verified=False,
             full_physical_residual_improvement_verified=False,
             full_PDE_corrections_verified=False, smooth_force_verified=False,
             complete_field_energy_verified=False,
             complete_field_blowup_lower_bound_verified=False, blowup_verified=False,
             formal_proof_assistant_certificate=False)
RESOLVED = ['conditional_actual_normalized_background_phase_and_moving_frame_comparison',
            'conditional_uniform_pulse_radial_lower_bound_and_all_fixed_derivative_bounds',
            'conditional_primary_covariance_columns_positive_weights_and_flat_edge_roots',
            'conditional_signed_weighted_inverse_at_fixed_primary_amplitudes',
            'conditional_full_linear_wave_remainder_gain_with_supported_sources',
            'conditional_real_primary_covariance_relative_gain_after_exact_curls']
REMAINING = ['physical_mean_velocity_pressure_and_integral_compatibility_corrections',
             'all_nonlinear_old_new_interactions_and_complete_residual_cycle_gain',
             'full_divergence_preserving_correction_sum_and_localization_initialization',
             'smooth_compact_forcing_through_T_complete_field_energy_and_blowup_lower_bound',
             'independent_review_of_inherited_continuum_bounds_and_global_auxiliary_geometry']
REVIEW_NAMES = ['phase_and_chart', 'band_and_frame', 'thresholds_and_barriers',
                'huge_band_pulse_fixture', 'cone_and_edges', 'derivatives_and_cutoff',
                'common_support_path', 'full_linear_operator', 'physical_partition_scaling']
GATE_NAMES = REVIEW_NAMES+['all_inherited_bytes_preserved', 'conditional_blocks_not_promoted_to_complete_solution']


def expected_protocol():
    return dict(status=STATUS, parent_commit='9c490d8732da619500584458f9bf71f31079a8aa',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['4.26-4.27', '5.42', '6.1-6.32', '7.2-7.40', '9.1-9.2'],
                actual_viscosity=1,
                actual_h='fixed inherited h in (0,1/100), not numerically selected here',
                actual_leading_radial_frequency='unchanged inherited N=1+floor(H^32)',
                diagnostic_h='1/200', diagnostic_kappa='1/100000', diagnostic_precision=120,
                diagnostic_first_band=51600, diagnostic_band_offsets=[0, 4000, 12000],
                diagnostic_rectangle_radius='.04', diagnostic_representative=['1.4', '.12', '.62'],
                diagnostic_background='manufactured implicit-chart profile plus epsilon-squared perturbations',
                pulse_solver='Radau on ratio and log(z_plus/P), at v/L',
                covariance_window_in_gaussian_coordinates=[-20, 20],
                covariance_quadrature_is_truncated_diagnostic=True,
                maximum_mp_identity_error='1e-85', maximum_binary64_operator_relative_error=2e-10,
                refined_pulse_column_relative_tolerance=2e-8,
                scope='Conditional actual-band phase/frame, positive pulse and all fixed derivatives, primary covariance/edge inverse and full linear gain. Numeric data are manufactured, actual constants and qstar are not evaluated; physical means, nonlinear complete-cycle gain, forcing, energy and blow-up remain open.')


def proof_ledger():
    return dict(normal_error='C*(S^-1+S*(S^-3+epsilon+epsilon^2+k^-1))',
                value_smallness='S^2*(epsilon+epsilon^2+k^-1)<=1',
                integer_rounding_bound='|p-ptilde|<=1/k; not uniformly 1/(2k)',
                single_qstar_depends_on_fixed_value_gaps_not_derivative_order_or_finite_harmonics=True,
                homogeneous_coefficient_bound='C_I*S^(2*|I|)*P',
                localized_primary_bound='C_I*sqrt(epsilon)*S^(3*|I|+1/4)*sqrt(zeta)*delta^-b_I*P',
                covariance_matrix_derivative_bound='|D^j H|<=C_j*S^(2j-1/2)',
                inverse_derivative_bound='|D^j H^-1|<=C_j*S^(2j+1/2)',
                positive_squared_amplitude_lower_bound='y_sigma>=c*sqrt(S)*zeta',
                signed_inverse_requires_angular_and_auxiliary_independent_target=True,
                full_actual_velocity_signed_quadratic_remainder_retained=True,
                real_primary_covariance_bound='epsilon^(2-2*kappa)*S^(3*|I|+6)*zeta*delta^-b_I',
                conservative_linear_gain='1/2-3*kappa',
                full_linear_operator_includes_pressure_viscosity_and_radial_fast_multiplier=True,
                global_linear_gain_requires_original_source_harmonic_smooth_temporal_extension=True,
                primary_covariance_does_not_cancel_all_mean_tensor_terms=True,
                full_nonlinear_residual_gain_supplied=False)


def provenance():
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit',
               'intermediate_decay_audit', 'pulse_moment_audit', 'pulse_stress_audit',
               'post_pulse_stress_audit', 'heat_exterior_audit', 'axis_matching_audit',
               'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds',
               'repair_state_bounds', 'physical_residual_budget', 'first_correction_construction',
               'higher_jet_budget', 'annular_first_correction', 'background_summation', 'physical_wave_interface')
    checks = {}
    for folder in folders:
        record = json.loads((HERE.parent/folder/'evidence.json').read_text())
        if 'source_sha256' in record:
            inherited, root = record['source_sha256'], HERE.parent/folder
        elif 'source_hashes' in record:
            inherited, root = record['source_hashes'], PROJECT
        else:
            inherited, root = record.get('provenance', {}), None
        if not inherited:
            raise ValueError('Missing inherited hash ledger: '+folder)
        for path, digest in inherited.items():
            routed = root if root is not None else (PROJECT if path.startswith('results/') or path == 'GOAL.md' else HERE.parent)
            checks[folder+':'+path] = hashlib.sha256((routed/path).read_bytes()).hexdigest() == digest
    return checks


def hashes():
    names = ('README.md', 'requirements.txt', 'protocol.json', 'bounds.py', 'charts.py',
             'diagnostics.py', 'audit.py', 'test_audit.py')
    paths = [HERE/name for name in names]+[PROJECT/'GOAL.md', HERE.parent/'background_summation/evidence.json',
                                          HERE.parent/'physical_wave_interface/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def encode(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, mp.mpf):
        if not mp.isfinite(value):
            raise ValueError('Nonfinite high-precision result')
        return mp.nstr(value, 80)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite binary64 result')
        return float(value)
    return value


def generate():
    mp.mp.dps = 120
    protocol = json.loads((HERE/'protocol.json').read_text())
    if protocol != expected_protocol() or fixture_first_band() != 51600:
        raise ValueError('Frozen diagnostic protocol or scope changed')
    reviews = dict(phase_and_chart=phase_identity_review(), band_and_frame=band_frame_review(),
                   thresholds_and_barriers=threshold_and_barrier_review(), huge_band_pulse_fixture=pulse_covariance_review(),
                   cone_and_edges=cone_and_edge_review(), derivatives_and_cutoff=derivative_cutoff_review(),
                   common_support_path=support_review(), full_linear_operator=linear_remainder_review(),
                   physical_partition_scaling=partition_scaling_review())
    inherited = provenance()
    gates = {name: all(reviews[name]['checks'].values()) for name in REVIEW_NAMES}
    gates['all_inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['conditional_blocks_not_promoted_to_complete_solution'] = not any(FLAGS[k] for k in
        ('actual_global_background_numerically_evaluated', 'actual_profile_constants_or_actual_qstar_numerically_selected',
         'global_auxiliary_color_centers_and_rectangles_materialized', 'actual_mean_correctors_constructed',
         'full_nonlinear_interaction_and_residual_cycle_verified', 'full_physical_stress_wave_construction_verified',
         'full_physical_residual_improvement_verified', 'full_PDE_corrections_verified', 'smooth_force_verified',
         'complete_field_energy_verified', 'complete_field_blowup_lower_bound_verified', 'blowup_verified'))
    return encode(dict(status=STATUS, protocol=protocol, flags=FLAGS, conditional_proof=proof_ledger(),
                       reviews=reviews, inherited_provenance=inherited, source_hashes=hashes(), gates=gates,
                       resolved_obligations=RESOLVED, remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS or record.get('conditional_proof') != proof_ledger():
        raise ValueError('Conditional pulse and linear estimates are not a full nonlinear PDE construction')
    if record.get('protocol') != expected_protocol():
        raise ValueError('Frozen diagnostic scope differs')
    if record.get('resolved_obligations') != RESOLVED or record.get('remaining_obligations') != REMAINING:
        raise ValueError('Mean, nonlinear, forcing or blow-up obligations changed')
    if set(record.get('gates', {})) != set(GATE_NAMES) or not all(record['gates'].values()):
        raise ValueError('A required band/pulse/linear/provenance gate failed or disappeared')


def compare_records(actual, frozen, location='root'):
    if isinstance(actual, dict):
        if not isinstance(frozen, dict) or set(actual) != set(frozen):
            raise ValueError('Frozen fields differ at '+location)
        for key in actual:
            compare_records(actual[key], frozen[key], location+'.'+key)
    elif isinstance(actual, list):
        if not isinstance(frozen, list) or len(actual) != len(frozen):
            raise ValueError('Frozen list differs at '+location)
        for i, (a, b) in enumerate(zip(actual, frozen)):
            compare_records(a, b, location+'['+str(i)+']')
    elif location.startswith('root.reviews.') and isinstance(actual, str) and isinstance(frozen, str):
        if '/' in actual or '/' in frozen:
            if actual != frozen:
                raise ValueError('Frozen rational differs at '+location)
            return
        try:
            a, b = mp.mpf(actual), mp.mpf(frozen)
        except ValueError:
            if actual != frozen:
                raise ValueError('Frozen label differs at '+location)
        else:
            if not mp.isfinite(a) or not mp.isfinite(b) or abs(a-b) > max(mp.mpf('1e-75'), mp.mpf('1e-20')*max(abs(a), abs(b))):
                raise ValueError('Frozen high-precision diagnostic differs at '+location)
    elif location.startswith('root.reviews.') and isinstance(actual, float):
        if isinstance(frozen, bool) or not isinstance(frozen, (float, int)) or not math.isclose(actual, frozen, rel_tol=1e-7, abs_tol=2e-10):
            raise ValueError('Frozen binary64 diagnostic differs at '+location)
    elif type(actual) is not type(frozen) or actual != frozen:
        raise ValueError('Frozen metadata differs at '+location)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify-record', type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit('Output already exists; choose a fresh path')
    record = generate()
    validate(record)
    if args.verify_record:
        frozen = json.loads(args.verify_record.read_text())
        validate(frozen)
        compare_records(record, frozen)
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    print(STATUS)
    print(str(len(record['gates']))+' gates pass; actual-band pulse/primary/linear bounds remain conditional.')
    print('The physical mean solve and complete nonlinear residual cycle remain open; no blow-up proof.')


if __name__ == '__main__':
    main()
