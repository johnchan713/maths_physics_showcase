#!/usr/bin/env python3
"""Reproduce the scoped extension and reject promotion to full PDE blow-up."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
import numpy as np

from extension import HERE, PROJECT, geometry_ledger, construction_ledger, residual_gain_ledger
from checks import (kinematic_review, global_matching_review, memory_source_review,
                    moment_identity_review, remainder_review)

STATUS = 'conditional-first-global-background-coefficient-and-relative-gain'
FLAGS = dict(conditional_actual_collar_continuation_supplied=True,
             conditional_actual_global_first_background_coefficient_extended=True,
             conditional_actual_five_moments_restored=True,
             conditional_signed_first_stress_support_supplied=True,
             conditional_first_augmented_background_exponent_gain=True,
             zero_U_does_not_imply_zero_leading_memory=True,
             zero_axis_pressure_and_same_leading_profile_preserved=True,
             relies_on_inherited_analytic_collar_and_leading_exact_moments=True,
             actual_high_frequency_profile_numerically_resolved=False,
             effective_global_first_coefficient_constants_evaluated=False,
             full_physical_stress_wave_construction_verified=False,
             full_physical_residual_budget_verified=False,
             absolute_physical_residual_decay_verified=False,
             all_order_background_induction_verified=False,
             full_PDE_corrections_verified=False,
             smooth_force_verified=False, complete_field_energy_verified=False,
             blowup_verified=False, independent_foundation_review_completed=False,
             formal_proof_assistant_certificate=False)
RESOLVED = ['conditional_actual_first_coefficient_continuation_into_preserved_collar',
            'actual_I3_geometry_and_nonzero_radial_memory_source_identified',
            'conditional_actual_global_first_coefficient_and_five_linear_moments',
            'conditional_compact_signed_first_tangential_stress',
            'conditional_first_augmented_background_relative_exponent_gain']
REMAINING = ['all_order_background_coefficient_induction_and_derivative_constants',
             'divergence_preserving_background_summation',
             'physical_wave_covariance_mean_corrections_and_all_interactions',
             'complete_physical_residual_budget_and_localization',
             'smooth_compact_forcing_through_T',
             'complete_field_blowup_lower_bound_and_energy',
             'independent_review_of_inherited_continuum_estimates']


def expected_protocol():
    return dict(status=STATUS, parent_commit='4a5a098862090231505d1976378574d6034ee9cb',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['Theorem 4.6', '5.2-5.6', 'Lemma 5.2', 'A.9', 'Corollary B.6'],
                fixed_N='1+floor(H^32)', viscosity=1,
                collar_endpoint='Xa*exp(tc/4)', positive_order_patch='I3',
                pressure_support_edge='Xv, not the end of I3',
                leading_patch_zero_memory_required=False,
                diagnostics_precision=80, global_matching_samples=12,
                matrix_quadrature=dict(order=48, panels=8), independent_quadrature=dict(order=48, panels=12),
                maximum_normalized_moment_error='1e-35',
                maximum_midplane_omega_pressure_moment='1e-35',
                minimum_off_midplane_pressure_omission='.01', minimum_fifth_moment_omission='.1',
                minimum_local_pressure_tail='1e-8', maximum_memory_source_error=1e-10,
                maximum_conservative_identity_error=1e-10, maximum_remainder_identity_error=2e-12,
                scope='Conditional analytic first background extension and signed stress support, with a relative augmented-background exponent gain only. Manufactured diagnostics are not the actual profile; full waves, all-order summation, smooth force, energy and blow-up remain open.')


def provenance():
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit',
               'intermediate_decay_audit', 'pulse_moment_audit', 'pulse_stress_audit',
               'post_pulse_stress_audit', 'heat_exterior_audit', 'axis_matching_audit',
               'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds',
               'repair_state_bounds', 'physical_residual_budget', 'first_correction_construction',
               'higher_jet_budget')
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
    paths = [HERE/name for name in ('README.md', 'requirements.txt', 'protocol.json', 'extension.py',
                                    'fixtures.py', 'checks.py', 'audit.py', 'test_audit.py')]
    paths += [PROJECT/'GOAL.md']
    for folder in ('higher_jet_budget', 'first_correction_construction', 'axis_core_attachment',
                   'stress_realization_audit', 'heat_exterior_audit', 'repair_state_bounds'):
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
    geometry, construction = geometry_ledger(), construction_ledger()
    reviews = dict(kinematics=kinematic_review(), global_matching=global_matching_review(),
                   memory_source=memory_source_review(), conservative_moments=moment_identity_review(),
                   finite_remainder=remainder_review())
    inherited = provenance()
    gains = [residual_gain_ledger(derivative_order=m) for m in (0, 1, 3)]
    gates = dict(exact_separated_support_geometry=all(geometry['checks'].values()),
                 construction_preserves_axis_pressure_and_source=all(construction['checks'].values()),
                 divergence_and_nonzero_memory=all(reviews['kinematics']['checks'].values()),
                 five_moments_and_pressure_omission_controls=all(reviews['global_matching']['checks'].values()),
                 independent_radial_memory_source=all(reviews['memory_source']['checks'].values()),
                 conservative_stress_tail_identities=all(reviews['conservative_moments']['checks'].values()),
                 quadratic_and_viscous_remainder=all(reviews['finite_remainder']['checks'].values()),
                 relative_gain_still_has_singular_physical_bound=all(x['relative_gain'] > 0 and x['sufficient_all_component_exponent'] < 0 for x in gains),
                 inherited_bytes_preserved=bool(inherited) and all(inherited.values()),
                 full_solution_not_promoted=not any(FLAGS[k] for k in
                     ('full_physical_stress_wave_construction_verified', 'full_physical_residual_budget_verified',
                      'absolute_physical_residual_decay_verified', 'full_PDE_corrections_verified',
                      'smooth_force_verified', 'complete_field_energy_verified', 'blowup_verified')))
    return encode(dict(status=STATUS, flags=FLAGS, protocol=protocol, geometry=geometry,
                       conditional_construction=construction, residual_gain=gains, reviews=reviews,
                       inherited_provenance=inherited, source_hashes=hashes(), gates=gates,
                       resolved_obligations=RESOLVED, remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS:
        raise ValueError('A conditional background extension is not a full PDE or blow-up proof')
    if record.get('resolved_obligations') != RESOLVED or record.get('remaining_obligations') != REMAINING:
        raise ValueError('The remaining wave, summation or forcing obligations changed')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A support, residual identity, omission control or provenance gate failed')


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
    result = generate()
    validate(result)
    if args.verify_record:
        frozen = json.loads(args.verify_record.read_text())
        validate(frozen)
        compare_records(result, frozen)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(STATUS)
    print(str(len(result['gates']))+' gates pass; first global background coefficient is conditional.')
    print('Relative augmented-background gain only: full waves, smooth force and blow-up remain open.')


if __name__ == '__main__':
    main()
