#!/usr/bin/env python3
"""Conditional background summation, with the retained physical-wave gap."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
import numpy as np

from hierarchy import HERE, PROJECT, induction_ledger
from summation import proof_ledger
from checks import (matrix_review, moment_review, cutoff_review, scale_review,
                    scalar_summation_review, exponent_review)

STATUS = 'conditional-all-order-background-with-flat-augmented-residual'
FLAGS = dict(conditional_all_order_background_induction_supplied=True,
             conditional_divergence_preserving_background_sum_supplied=True,
             conditional_augmented_background_error_flat_in_all_fixed_derivatives=True,
             conditional_background_swirl_lower_bound_survives_all_background_orders=True,
             zero_axis_pressure_mean_patch_velocity_and_heat_exterior_preserved=True,
             shifted_radial_axial_viscosity_source_retained=True,
             relies_on_inherited_leading_profile_collar_moments_and_heat_bounds=True,
             actual_high_frequency_profile_numerically_resolved=False,
             actual_coefficient_sequence_or_cutoff_scales_numerically_evaluated=False,
             uncut_formal_power_series_convergence_asserted=False,
             uniform_corrected_stress_cone_verified=False,
             bare_background_momentum_residual_flat=False,
             full_physical_stress_wave_construction_verified=False,
             full_physical_residual_budget_verified=False,
             full_PDE_corrections_verified=False,
             smooth_force_verified=False, complete_field_energy_verified=False,
             complete_field_blowup_lower_bound_verified=False,
             blowup_verified=False, independent_foundation_review_completed=False,
             formal_proof_assistant_certificate=False)
RESOLVED = ['conditional_general_order_regular_background_system_and_shifted_sources',
            'conditional_common_collar_induction_and_five_linear_moment_repairs',
            'conditional_common_higher_stress_support_and_fixed_derivative_losses',
            'conditional_locally_finite_divergence_preserving_asymptotic_background_sum',
            'conditional_flat_augmented_background_residual',
            'conditional_summed_background_core_swirl_lower_bound']
REMAINING = ['full_physical_wave_covariance_and_signed_stress_realization',
             'physical_mean_corrections_and_all_cross_interactions',
             'full_field_residual_improvement_and_divergence_preserving_summation',
             'localization_initialization_and_smooth_compact_forcing_through_T',
             'complete_field_blowup_lower_bound_and_energy',
             'independent_review_of_inherited_continuum_estimates']


def expected_protocol():
    return dict(status=STATUS, parent_commit='305b9288a86877ba994514e514835a0c8e0f036d',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['5.1-5.6', 'Lemmas 5.1-5.4', '5.41-5.46'],
                fixed_N='1+floor(H^32)', actual_viscosity=1,
                same_collar='0<=X<=Xa*exp(tc/4)', numerical_general_orders=[1, 2, 3, 5],
                matrix_samples=216, moment_samples=48, cutoff_samples=96,
                maximum_system_relative_error=3e-11, maximum_first_system_error=3e-12,
                maximum_moment_identity_error=3e-11, maximum_cutoff_divergence_error=3e-11,
                maximum_curl_value_gradient_error=3e-11, minimum_omission_gap=1e-5,
                scale_check_precisions=[80, 110], scalar_diagnostic_precision=80,
                scalar_cutoff_active_indices=[4, 8, 12, 16],
                scalar_ODE='q^2*f_prime-f+q+f^2',
                finite_residual_loss='3/2+2h+m, independent of truncation',
                scope='Conditional all-order background, common support, divergence-preserving asymptotic sum and flat augmented error. The retained stress is not a full physical wave covariance; the bare residual, full field, smooth forcing, energy and blow-up remain open.')


def provenance():
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit',
               'intermediate_decay_audit', 'pulse_moment_audit', 'pulse_stress_audit',
               'post_pulse_stress_audit', 'heat_exterior_audit', 'axis_matching_audit',
               'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds',
               'repair_state_bounds', 'physical_residual_budget', 'first_correction_construction',
               'higher_jet_budget', 'annular_first_correction')
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
    paths = [HERE/name for name in ('README.md', 'requirements.txt', 'protocol.json', 'hierarchy.py',
                                    'summation.py', 'checks.py', 'audit.py', 'test_audit.py')]
    paths += [PROJECT/'GOAL.md']
    for folder in ('annular_first_correction', 'higher_jet_budget', 'first_correction_construction',
                   'axis_core_attachment', 'stress_realization_audit', 'heat_exterior_audit', 'physical_residual_budget'):
        paths.append(HERE.parent/folder/'evidence.json')
    return {str(path.relative_to(PROJECT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def encode(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, mp.mpf):
        if not mp.isfinite(value):
            raise ValueError('Nonfinite multiprecision result')
        return mp.nstr(value, 80)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite numerical diagnostic')
        return float(value)
    return value


def generate():
    mp.mp.dps = 80
    protocol = json.loads((HERE/'protocol.json').read_text())
    if protocol != expected_protocol():
        raise ValueError('Frozen construction, diagnostic thresholds or scope changed')
    reviews = dict(matrices=matrix_review(), conservative_moments=moment_review(),
                   potential_cutoffs=cutoff_review(), logarithmic_scales=scale_review(),
                   scalar_summation=scalar_summation_review(), exponents=exponent_review())
    inherited = provenance()
    gates = {name: all(record['checks'].values()) for name, record in reviews.items()}
    gates['inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['conditional_background_not_promoted_to_full_solution'] = not any(FLAGS[k] for k in
        ('bare_background_momentum_residual_flat', 'uniform_corrected_stress_cone_verified',
         'full_physical_stress_wave_construction_verified', 'full_physical_residual_budget_verified',
         'full_PDE_corrections_verified', 'smooth_force_verified', 'complete_field_energy_verified',
         'complete_field_blowup_lower_bound_verified', 'blowup_verified'))
    return encode(dict(status=STATUS, flags=FLAGS, protocol=protocol,
                       conditional_induction=induction_ledger(), conditional_summation=proof_ledger(),
                       reviews=reviews, inherited_provenance=inherited, source_hashes=hashes(), gates=gates,
                       resolved_obligations=RESOLVED, remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS:
        raise ValueError('Flat augmented background error is not a full PDE or blow-up proof')
    if record.get('resolved_obligations') != RESOLVED or record.get('remaining_obligations') != REMAINING:
        raise ValueError('The retained stress, wave or forcing obligations changed')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A hierarchy, curl, summation or provenance gate failed')


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
        # Exact rational ledgers are metadata, not approximate decimals.
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
            if not mp.isfinite(a) or not mp.isfinite(b) or abs(a-b) > max(mp.mpf('1e-60'), mp.mpf('1e-12')*max(abs(a), abs(b))):
                raise ValueError('Frozen multiprecision diagnostic differs at '+location)
    elif location.startswith('root.reviews.') and isinstance(actual, float):
        if isinstance(frozen, bool) or not isinstance(frozen, (float, int)) or not math.isclose(actual, frozen, rel_tol=1e-8, abs_tol=2e-13):
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
    print(str(len(record['gates']))+' gates pass; all-order background and summation are conditional.')
    print('The retained stress remains: bare residual, physical waves, force, energy and blow-up stay open.')


if __name__ == '__main__':
    main()
