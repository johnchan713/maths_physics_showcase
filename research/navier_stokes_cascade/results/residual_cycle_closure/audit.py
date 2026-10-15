#!/usr/bin/env python3
"""Conditional initial state and complete finite-cycle order gain; no blow-up."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from cycle_bounds import (Affine, KAPPA, MIN_B, STEP, cycle_certificate,
                          initialization, physical_loss, physical_residual_power)
from interaction_checks import (covariance_increment_review, full_increment_review,
                                longitudinal_and_transport_review)
from support_checks import sharp_loss_controls, support_and_flat_review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'conditional-initial-state-and-complete-finite-cycle-order-gain'
FLAGS = dict(conditional_initial_corrected_state_supplied=True,
             complete_curl_longitudinal_and_nonlinear_estimates_derived=True,
             source_closure_proved_given_inherited_support_geometry=True,
             conditional_complete_finite_cycle_order_gain_supplied=True,
             conditional_common_domain_and_stage_independent_physical_jet_loss_supplied=True,
             actual_global_profile_and_field_numerically_evaluated=False,
             inherited_foundation_and_global_auxiliary_geometry_independently_verified=False,
             complete_infinite_wave_and_mean_potential_sum_verified=False,
             unconditional_complete_residual_construction_verified=False,
             low_stage_physical_residual_proved_to_decay=False,
             fixed_scale_residual_norm_contraction_verified=False,
             admissible_smooth_compact_force_verified=False,
             complete_field_energy_verified=False,
             complete_field_blowup_lower_bound_verified=False, blowup_verified=False,
             diagnostic_success_is_not_a_continuum_proof=True,
             formal_proof_assistant_certificate=False)
REVIEW_NAMES = ('initial_state_orders', 'all_stage_cycle_orders', 'longitudinal_and_transport',
                'complete_signed_covariance', 'full_residual_increment',
                'source_support_and_flat_memory', 'sharp_missing_hypothesis_controls',
                'fixed_physical_jet_loss')
GATE_NAMES = REVIEW_NAMES+('all_inherited_bytes_preserved', 'conditional_finite_result_not_promoted')
REMAINING = ('full_wave_and_mean_potential_cutoff_sum_with_all_curl_commutators',
             'smooth_realized_residual_through_singular_time_and_localization_initialization',
             'complete_field_energy_and_surviving_velocity_lower_bound',
             'independent_inherited_profile_background_and_global_auxiliary_geometry_review')


def expected_protocol():
    return dict(status=STATUS, parent_commit='73116848b1e53cad36e3c0e314aaaac452518a2e',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['6.8-6.32', '7.39-7.42', '9.1-9.10', 'Propositions 9.3-9.6'],
                actual_profile_parameters='unchanged inherited N, h, positive lambda and fixed inverse maps',
                fixed_kappa='1/100000', minimum_B='7/10', cycle_order_gain='1/10',
                primary_real_covariance_remainder_order='2-2*kappa',
                total_wave_remainder_order='17/25', total_mean_order='9/10',
                all_stage_inequalities='affine half-line proof on B>=7/10; not a stage scan',
                diagnostic_frequencies=[16, 64, 256, 1024],
                diagnostic_covariance_frequencies=[11, 29],
                diagnostic_jet_points=[[1.12, -.4, -.19, .36], [1.3, .27, .2, .47], [1.66, -.9, .31, .64]],
                diagnostic_physical_h='1/200', diagnostic_float_relative_tolerance=2e-10,
                diagnostic_data='manufactured compact curls and means plus separate Cartesian missing-hypothesis examples',
                scope='Conditional initialization, source closure and complete finite-cycle order gain with a fixed physical jet loss. Actual global evaluation, independent foundation review, infinite realization, forcing, energy and blow-up remain open.')


def proof_ledger():
    return dict(longitudinal_identity='n.a=i*div_star(a)/(k*m) for the COMPLETE curl amplitude',
                same_label_wave_wave_order='alpha+alpha_prime-kappa',
                wave_mean_order='alpha+mu-1/2; radial mean has order mu+1',
                distinct_closed_label_supports_are_an_inherited_premise=True,
                mean_operations_may_fill_auxiliary_torus=True,
                subsequent_nonzero_sources_retain_a_wave_factor_and_its_P_envelope=True,
                angular_harmonic_bound='M_next<=2*max(M,1); finite at every fixed stage',
                pressure_recomputed_after_every_velocity_change=True,
                moving_axial_pressure_product_retained='epsilon*partial_Z(Jz+c_rho*P)',
                old_new_and_new_square_covariances_retained=True,
                flat_velocities_remain_in_all_actual_products=True,
                initial_orders=['1-3*kappa', '2-3*kappa', '1.9-3*kappa'],
                stage_orders=['B=7/10+j/10', 'C=B+1/2'],
                complete_cycle_gain='1/10 in declared normalized residual orders',
                all_finite_states_use_fixed_primary_denominators_and_one_domain=True,
                physical_jet_loss='K_m=2*A+3/2+2*m; independent of finite stage',
                physical_bound='C_jm*q^(h*B_j-K_m); j fixed before q tends to zero',
                no_uniform_stage_constant_or_raw_series_convergence_claim=True)


def inherited_provenance():
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit',
               'intermediate_decay_audit', 'pulse_moment_audit', 'pulse_stress_audit',
               'post_pulse_stress_audit', 'heat_exterior_audit', 'axis_matching_audit',
               'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds',
               'repair_state_bounds', 'physical_residual_budget', 'first_correction_construction',
               'higher_jet_budget', 'annular_first_correction', 'background_summation',
               'physical_wave_interface', 'actual_band_pulse_bounds', 'compact_mean_corrections')
    checks = {}
    for folder in folders:
        record = json.loads((HERE.parent/folder/'evidence.json').read_text())
        if 'source_sha256' in record:
            hashes, root = record['source_sha256'], HERE.parent/folder
        elif 'source_hashes' in record:
            hashes, root = record['source_hashes'], PROJECT
        else:
            hashes, root = record.get('provenance', {}), None
        if not hashes:
            raise ValueError('Missing inherited ledger: '+folder)
        for path, digest in hashes.items():
            routed = root if root is not None else (PROJECT if path.startswith('results/') or path == 'GOAL.md' else HERE.parent)
            checks[folder+':'+path] = hashlib.sha256((routed/path).read_bytes()).hexdigest() == digest
    return checks


def source_hashes():
    names = ('README.md', 'requirements.txt', 'protocol.json', 'cycle_bounds.py',
             'interaction_checks.py', 'support_checks.py', 'audit.py', 'test_audit.py')
    paths = [HERE/name for name in names]+[PROJECT/'GOAL.md',
        HERE.parent/'compact_mean_corrections/evidence.json',
        HERE.parent/'actual_band_pulse_bounds/evidence.json',
        HERE.parent/'background_summation/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def encode(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, Affine):
        return dict(slope=str(value.slope), constant=str(value.constant))
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite diagnostic')
        return float(value)
    return value


def physical_review():
    h = Fraction(1, 200)
    rows = [dict(derivative_order=m, fixed_loss=physical_loss(h, m),
                 initial_power=physical_residual_power(0, m, h),
                 one_cycle_power=physical_residual_power(1, m, h)) for m in (0, 1, 3, 8)]
    return dict(rows=rows, phase_time_cost_upper_bound_at_h_limit=1+Fraction(5, 2)*Fraction(1, 100),
                power_slope_in_stage=h*STEP,
                checks=dict(physical_derivative_cost_is_below_two=1+Fraction(5, 2)*Fraction(1, 100) < 2,
                            stage_independent_jet_loss=all(r['one_cycle_power']-r['initial_power'] == h*STEP for r in rows),
                            normalized_gain_not_promoted_to_low_stage_absolute_decay=rows[0]['initial_power'] < 0,
                            strict_order_gain_is_positive=h*STEP > 0))


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != expected_protocol():
        raise ValueError('Frozen scope or parameters changed')
    initial, cycle = initialization(), cycle_certificate()
    initial['checks'] = dict(all_initial_state_margins_strict=min(initial['margins'].values()) > 0)
    cycle['checks'] = dict(all_half_line_margins_nonnegative=min(cycle['all_half_line_margins'].values()) >= 0,
                          declared_gain_is_positive=cycle['requested_gain'] > 0)
    reviews = dict(initial_state_orders=initial, all_stage_cycle_orders=cycle,
                   longitudinal_and_transport=longitudinal_and_transport_review(),
                   complete_signed_covariance=covariance_increment_review(),
                   full_residual_increment=full_increment_review(),
                   source_support_and_flat_memory=support_and_flat_review(),
                   sharp_missing_hypothesis_controls=sharp_loss_controls(),
                   fixed_physical_jet_loss=physical_review())
    inherited = inherited_provenance()
    gates = {name: all(reviews[name]['checks'].values()) for name in REVIEW_NAMES}
    gates['all_inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['conditional_finite_result_not_promoted'] = not any(FLAGS[k] for k in
        ('actual_global_profile_and_field_numerically_evaluated',
         'inherited_foundation_and_global_auxiliary_geometry_independently_verified',
         'complete_infinite_wave_and_mean_potential_sum_verified',
         'unconditional_complete_residual_construction_verified',
         'low_stage_physical_residual_proved_to_decay', 'fixed_scale_residual_norm_contraction_verified',
         'admissible_smooth_compact_force_verified', 'complete_field_energy_verified',
         'complete_field_blowup_lower_bound_verified', 'blowup_verified'))
    return encode(dict(status=STATUS, protocol=expected_protocol(), flags=FLAGS,
                       conditional_proof=proof_ledger(), reviews=reviews,
                       inherited_provenance=inherited, source_hashes=source_hashes(),
                       gates=gates, remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS or record.get('conditional_proof') != proof_ledger():
        raise ValueError('The conditional finite-cycle theorem is not an unconditional blow-up construction')
    if record.get('protocol') != expected_protocol() or record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('Frozen scope or remaining obligations changed')
    if set(record.get('gates', {})) != set(GATE_NAMES) or not all(record['gates'].values()):
        raise ValueError('A required cycle, interaction, support or provenance gate failed')


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
    elif location.startswith('root.reviews.') and isinstance(actual, float):
        if isinstance(frozen, bool) or not isinstance(frozen, (float, int)) or not math.isclose(actual, frozen, rel_tol=1e-7, abs_tol=2e-10):
            raise ValueError('Frozen diagnostic differs at '+location)
    elif type(actual) is not type(frozen) or actual != frozen:
        raise ValueError('Frozen exact arithmetic or metadata differs at '+location)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
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
    args.output.write_text(json.dumps(record, sort_keys=True, indent=2)+'\n')
    print(STATUS)
    print(str(len(record['gates']))+' gates pass; the complete finite-cycle gain remains conditional.')
    print('Infinite realization, independent foundation review, smooth forcing, energy and blow-up remain open.')


if __name__ == '__main__':
    main()
