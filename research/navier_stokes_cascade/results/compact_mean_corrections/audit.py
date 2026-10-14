#!/usr/bin/env python3
"""Reproduce compact mean identities and a conditional compatibility gain."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
import numpy as np
from reviews import (defect_increment_review, integrated_mean_review, moment_map_review,
                     physical_mean_review, torus_and_radial_review)

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'conditional-compact-mean-inverses-and-compatibility-defect-gain'
FLAGS = dict(complete_six_entry_mean_balance_and_axial_pressure_moment_derived=True,
             conditional_compact_phase_path_radial_inverse_supplied=True,
             conditional_pressure_and_divergence_free_mean_maps_supplied=True,
             conditional_zero_auxiliary_mean_temporal_inverse_supplied=True,
             conditional_fixed_signed_five_equation_mean_inverse_supplied=True,
             conditional_slow_mean_update_compatibility_defect_gain_supplied=True,
             actual_global_background_and_mean_fields_numerically_evaluated=False,
             actual_primary_mean_initialization_verified=False,
             all_nonlinear_wave_mean_interactions_and_source_closure_verified=False,
             full_residual_cycle_gain_verified=False,
             smooth_force_verified=False, complete_field_energy_verified=False,
             complete_field_blowup_lower_bound_verified=False, blowup_verified=False,
             inherited_continuum_foundation_and_support_geometry_remain_conditional=True,
             independent_foundation_review_completed=False,
             finite_quadrature_and_stencil_diagnostics_are_interval_certificates=False,
             formal_proof_assistant_certificate=False)
REVIEW_NAMES = ('full_mean_operator', 'compact_radial_and_temporal_maps',
                'pressure_and_integrated_constraints', 'five_equation_inverse',
                'same_field_nonlinear_defect_changes')
GATE_NAMES = REVIEW_NAMES+('all_inherited_source_bytes_preserved', 'complete_solution_not_inferred')
REMAINING = ('actual_primary_mean_initialization_and_complete_nonlinear_cycle_gain',
             'nonzero_harmonic_source_support_and_gaussian_envelopes_after_all_mean_operations',
             'uniform_fixed_physical_derivative_loss_and_full_potential_sum',
             'smooth_compact_forcing_complete_field_energy_and_blowup_lower_bound',
             'independent_inherited_foundation_and_global_auxiliary_geometry_review')


def expected_protocol():
    return dict(status=STATUS, parent_commit='6498412c39cf120bea94b90c1f25fc33bf2cea8a',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['6.2', '6.6', '8.1-8.27'],
                actual_profile_frequency_and_parameters='unchanged inherited N, h and positive lambda',
                diagnostic_precision=80, diagnostic_kappa='1/100000',
                diagnostic_moment_lambdas=['.2', '.05', '.001'],
                diagnostic_patch_amplitudes=['1.3', '-.8'],
                diagnostic_radial_inverse_exponent=2,
                diagnostic_radial_source='degree-16 polynomial in U with C7 zero extension; not the C-infinity continuum source',
                diagnostic_radial_scales=[1, 8, 64, 512],
                diagnostic_normalized_epsilon=.13,
                diagnostic_pressure_nodes_per_support_segment=96,
                diagnostic_pressure_z_five_point_step=.002,
                diagnostic_refined_pressure_nodes_per_support_segment=128,
                diagnostic_refined_pressure_z_five_point_step=.001,
                physical_identity_relative_tolerance=2e-10,
                integrated_identity_absolute_tolerance=1e-8,
                diagnostic_data='manufactured compact curls and means, finite Fourier modes and fixed bump moments',
                scope='Conditional compact mean maps and slow compatibility-defect gain. Actual global mean initialization, complete residual cycle, forcing, energy and blow-up remain open.')


def proof_ledger():
    return dict(radial_graph_inverse='D_e T_e f=f-A_e f, e=0,1,2; phase path retained',
                radial_inverse_preserves_input_epsilon_order=True,
                radial_cutoff_flatness_requires_zero_weighted_auxiliary_mean_and_all_fixed_input_derivatives=True,
                radial_cutoff_remainder_order='alpha+p*kappa for every fixed integer p>=1',
                temporal_inverse_requires_zero_auxiliary_mean_at_each_slow_radial_point=True,
                temporal_inverse_derivative_loss=4,
                mean_maps_can_fill_auxiliary_torus=True,
                pressure_radial_residual='-rho*P-A_0(gr-rho*P)',
                complete_axial_integrated_flux='epsilon*partial_Z(Jz+c_rho*P)',
                moving_c_rho_derivative_is_retained=True,
                five_map_uses_fixed_positive_lambda_and_nonzero_patch_amplitude=True,
                five_map_has_no_sign_or_small_target_restriction=True,
                five_map_does_not_set_nonlinear_defects_to_zero=True,
                assumed_current_mean_orders=['0.9', '0.9', '1.9'],
                assumed_target_order='alpha>=0.9',
                induced_radial_increment_order='alpha+1',
                conditional_updated_compatibility_defect_order='alpha+0.9-2*kappa',
                full_old_new_wave_and_mean_cycle_gain_supplied=False)


def inherited_provenance():
    # Preserve the exact routing used in the saved parent audit. Read each
    # expected digest from its frozen source ledger; do not trust pass flags.
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit',
               'intermediate_decay_audit', 'pulse_moment_audit', 'pulse_stress_audit',
               'post_pulse_stress_audit', 'heat_exterior_audit', 'axis_matching_audit',
               'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds',
               'repair_state_bounds', 'physical_residual_budget', 'first_correction_construction',
               'higher_jet_budget', 'annular_first_correction', 'background_summation',
               'physical_wave_interface', 'actual_band_pulse_bounds')
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
            raise ValueError('Missing inherited source hashes: '+folder)
        for path, digest in hashes.items():
            routed = root if root is not None else (PROJECT if path.startswith('results/') or path == 'GOAL.md' else HERE.parent)
            checks[folder+':'+path] = hashlib.sha256((routed/path).read_bytes()).hexdigest() == digest
    return checks


def source_hashes():
    names = ('README.md', 'requirements.txt', 'protocol.json', 'operators.py',
             'mean_system.py', 'reviews.py', 'audit.py', 'test_audit.py')
    paths = [HERE/name for name in names]+[PROJECT/'GOAL.md',
             HERE.parent/'actual_band_pulse_bounds/evidence.json',
             HERE.parent/'physical_wave_interface/evidence.json',
             HERE.parent/'background_summation/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def encode(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, mp.mpf):
        if not mp.isfinite(value):
            raise ValueError('Nonfinite high-precision result')
        return mp.nstr(value, 70)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite binary64 result')
        return float(value)
    return value


def generate():
    mp.mp.dps = 80
    if json.loads((HERE/'protocol.json').read_text()) != expected_protocol():
        raise ValueError('Frozen protocol or mathematical scope changed')
    reviews = dict(full_mean_operator=physical_mean_review(),
                   compact_radial_and_temporal_maps=torus_and_radial_review(),
                   pressure_and_integrated_constraints=integrated_mean_review(),
                   five_equation_inverse=moment_map_review(),
                   same_field_nonlinear_defect_changes=defect_increment_review())
    provenance = inherited_provenance()
    gates = {name: all(reviews[name]['checks'].values()) for name in REVIEW_NAMES}
    gates['all_inherited_source_bytes_preserved'] = bool(provenance) and all(provenance.values())
    gates['complete_solution_not_inferred'] = not any(FLAGS[k] for k in
        ('actual_global_background_and_mean_fields_numerically_evaluated',
         'actual_primary_mean_initialization_verified',
         'all_nonlinear_wave_mean_interactions_and_source_closure_verified',
         'full_residual_cycle_gain_verified', 'smooth_force_verified', 'complete_field_energy_verified',
         'complete_field_blowup_lower_bound_verified', 'blowup_verified'))
    return encode(dict(status=STATUS, protocol=expected_protocol(), flags=FLAGS,
                       conditional_proof=proof_ledger(), reviews=reviews,
                       inherited_provenance=provenance, source_hashes=source_hashes(),
                       gates=gates, remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS or record.get('conditional_proof') != proof_ledger():
        raise ValueError('A compact mean estimate must not be promoted to a complete residual or blow-up proof')
    if record.get('protocol') != expected_protocol() or record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('Frozen protocol or open obligations changed')
    if set(record.get('gates', {})) != set(GATE_NAMES) or not all(record['gates'].values()):
        raise ValueError('A required compact mean, increment or provenance gate failed')


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
            if not mp.isfinite(a) or not mp.isfinite(b) or abs(a-b) > max(mp.mpf('1e-50'), mp.mpf('1e-35')*max(abs(a), abs(b))):
                raise ValueError('Frozen high-precision diagnostic differs at '+location)
    elif location.startswith('root.reviews.') and isinstance(actual, float):
        if isinstance(frozen, bool) or not isinstance(frozen, (float, int)) or not math.isclose(actual, frozen, rel_tol=1e-7, abs_tol=2e-10):
            raise ValueError('Frozen binary64 diagnostic differs at '+location)
    elif type(actual) is not type(frozen) or actual != frozen:
        raise ValueError('Frozen metadata differs at '+location)


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
    print(str(len(record['gates']))+' gates pass; compact mean and compatibility-defect gains remain conditional.')
    print('Actual initialization, complete nonlinear residual gain, forcing, energy and blow-up remain open.')


if __name__ == '__main__':
    main()
