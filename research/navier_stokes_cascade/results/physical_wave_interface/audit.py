#!/usr/bin/env python3
"""Physical curl/covariance/residual interface with actual pulse bounds open."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from checks import (POINTS, FREQUENCIES, coordinate_review, covariance_review,
                    residual_review, mean_pressure_review)
from pulses import pulse_review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'physical-wave-interface-with-open-actual-pulse-and-mean-bounds'
FLAGS = dict(exact_real_curl_and_full_covariance_identity_supplied=True,
             full_cylindrical_mean_residual_interface_supplied=True,
             conditional_moving_normal_pulse_propagator_argument_supplied=True,
             compact_manufactured_physical_wave_checked_in_two_coordinates=True,
             normalized_reference_shear_viscosity_pulse_checked=True,
             generic_kinematic_wave_residual_obstruction_supplied=True,
             actual_band_phase_and_pulse_uniform_bounds_verified=False,
             actual_band_support_descent_and_all_cross_interactions_verified=False,
             uniform_actual_signed_inverse_and_edge_weights_verified=False,
             actual_mean_correctors_constructed=False,
             full_physical_stress_wave_construction_verified=False,
             full_physical_residual_improvement_verified=False,
             full_PDE_corrections_verified=False, smooth_force_verified=False,
             complete_field_energy_verified=False,
             complete_field_blowup_lower_bound_verified=False, blowup_verified=False,
             independent_foundation_review_completed=False,
             formal_proof_assistant_certificate=False)
RESOLVED = ['exact_real_curl_carrier_and_quadratic_full_covariance_correction',
            'exact_full_mean_tensor_divergence_and_radial_pressure_defect',
            'exact_signed_increment_identity_retaining_all_remainders',
            'conditional_moving_normal_principal_equation_and_propagator_lemma',
            'manufactured_physical_residual_and_reference_pulse_crosschecks',
            'kinematic_matching_alone_has_no_full_residual_improvement']
REMAINING = ['uniform_actual_band_phase_frames_pulses_and_all_fixed_derivatives',
             'actual_positive_primary_and_signed_higher_stress_with_edge_weights',
             'actual_mean_corrections_integral_constraints_and_all_interactions',
             'full_physical_residual_improvement_and_divergence_preserving_sum',
             'localization_initialization_and_smooth_compact_forcing_through_T',
             'complete_field_energy_and_blowup_lower_bound',
             'independent_review_of_inherited_continuum_estimates']


def expected_protocol():
    return dict(status=STATUS, parent_commit='25f59628da5bc7c2efb66816075ddd98d7e340c1',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['7.2-7.16', '7.24-7.41', '8.3', '9.1-9.2'],
                coordinates=['r', 'theta', 'z', 't'], potential_Taylor_order=3,
                actual_viscosity=1, manufactured_points=[list(p) for p in POINTS],
                covariance_frequencies=FREQUENCIES, physical_coordinate_frequencies=[8, 23],
                fixture_angular_multipliers=[1, 2], angular_quadrature_nodes='8*k+1',
                reference_pulse_lengths=[40, 100, 240], reference_epsilon_times_k_squared=1,
                maximum_coordinate_covariance_residual_relative_error=2e-10,
                maximum_reference_constraint_principal_energy_error=2e-9,
                scope='Exact curl, covariance, mean/residual and conditional pulse identities; compact manufactured physical waves and reference ODE diagnostics. Actual band-uniform pulse bounds, mean correctors, full residual improvement, forcing, energy and blow-up remain open.')


def provenance():
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit',
               'intermediate_decay_audit', 'pulse_moment_audit', 'pulse_stress_audit',
               'post_pulse_stress_audit', 'heat_exterior_audit', 'axis_matching_audit',
               'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds',
               'repair_state_bounds', 'physical_residual_budget', 'first_correction_construction',
               'higher_jet_budget', 'annular_first_correction', 'background_summation')
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
    names = ('README.md', 'requirements.txt', 'protocol.json', 'jets.py', 'waves.py',
             'pulses.py', 'checks.py', 'audit.py', 'test_audit.py')
    paths = [HERE/name for name in names]+[PROJECT/'GOAL.md', HERE.parent/'background_summation/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def encode(value):
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite wave diagnostic')
        return float(value)
    return value


def generate():
    protocol = json.loads((HERE/'protocol.json').read_text())
    if protocol != expected_protocol():
        raise ValueError('Frozen wave fixture, thresholds or scope changed')
    reviews = dict(physical_coordinates=coordinate_review(), full_covariance=covariance_review(),
                   full_residual=residual_review(), mean_pressure=mean_pressure_review(),
                   reference_pulses=pulse_review())
    inherited = provenance()
    gates = {name: all(record['checks'].values()) for name, record in reviews.items()}
    gates['all_inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['reference_and_identity_results_not_promoted_to_actual_waves'] = not any(
        FLAGS[k] for k in FLAGS if k.endswith('_verified') or k == 'actual_mean_correctors_constructed')
    return encode(dict(status=STATUS, protocol=protocol, flags=FLAGS, reviews=reviews,
                       inherited_provenance=inherited, source_hashes=hashes(), gates=gates,
                       resolved_obligations=RESOLVED, remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS:
        raise ValueError('Interface identities and reference diagnostics are not an actual PDE construction')
    if record.get('protocol') != expected_protocol():
        raise ValueError('Frozen wave protocol differs')
    if record.get('resolved_obligations') != RESOLVED or record.get('remaining_obligations') != REMAINING:
        raise ValueError('Actual pulse, mean, residual or blow-up obligations changed')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A wave identity, diagnostic or provenance gate failed')


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
            raise ValueError('Frozen numerical diagnostic differs at '+location)
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
    print(str(len(record['gates']))+' gates pass; actual band pulses and mean corrections remain open.')
    print('Covariance matching alone gives no full residual improvement or blow-up proof.')


if __name__ == '__main__':
    main()
