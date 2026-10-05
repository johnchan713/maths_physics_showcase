#!/usr/bin/env python3
"""Reproduce the conditional inner construction without promoting global blow-up."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
import numpy as np

from transfer import HERE, PROJECT, angular_transfer, inner_operator_certificate, norm_ledger, third_implicit_terms, lower, upper
from majorant import proof_ledger
from checks import matrix_review, sparse_review, recurrence_review, first_moment_review, implicit_review

STATUS = 'conditional-first-inner-background-correction-with-open-global-extension'
FLAGS = dict(actual_unmodulated_C3_transfer_supplied=True,
             actual_core_operator_envelope_supplied=True,
             conditional_actual_inner_background_coefficient_constructed=True,
             positive_order_five_moment_linear_map_derived=True,
             conditional_one_correction_core_swirl_lower_bound=True,
             relies_on_inherited_analytic_core=True,
             actual_post_modulation_C3_bound_verified=False,
             actual_global_first_coefficient_extension_verified=False,
             actual_annular_support_hypotheses_fully_verified=False,
             actual_high_frequency_profile_numerically_resolved=False,
             full_physical_residual_budget_verified=False,
             full_physical_stress_wave_construction_verified=False,
             full_PDE_corrections_verified=False,
             smooth_force_verified=False,
             complete_field_energy_verified=False,
             blowup_verified=False,
             independent_foundation_review_completed=False,
             formal_proof_assistant_certificate=False)
RESOLVED = ['conditional_actual_unmodulated_incoming_C3_transfer',
            'regular_six_variable_first_background_system',
            'conditional_actual_inner_operator_bound_and_Picard_construction',
            'conditional_positive_order_linear_five_moment_inverse',
            'conditional_one_inner_background_correction_swirl_lower_bound']
REMAINING = ['original_moment_C4_and_loop_C3_input_bounds',
             'actual_post_modulation_radial_and_C3_angular_bounds',
             'actual_first_coefficient_annular_extension_and_five_moment_support',
             'full_wave_covariance_and_all_nonlinear_interactions',
             'all_order_PDE_correction_summation',
             'localization_and_smooth_forcing_through_T',
             'complete_field_blowup_lower_bound_and_energy',
             'independent_review_of_inherited_continuum_estimates']


def expected_protocol():
    return dict(status=STATUS, parent_commit='5a94f20dfc1de702a507cc9764bfc4a68ff01ab8',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['5.1-5.6', 'Lemma 5.1', 'B.4', 'B.26', 'B.34-B.40'],
                fixed_N='1+floor(H^32)', actual_inner_interval='0<=X<=3/Lambda=(3/4)*Xa',
                actual_inner_viscosity=1, actual_inner_operator_envelope='C^16',
                actual_profile_numerically_resolved=False, interval_digits=[80, 110],
                matrix_sample_count=54, maximum_matrix_relative_error=2e-12,
                minimum_negative_control_gap=1e-5, manufactured_radial_degrees=[3, 5, 7],
                maximum_final_manufactured_residual=2e-11, minimum_manufactured_refinement=1000000,
                maximum_axis_slope_error=2e-12, positive_order_patch='I_3',
                patch_diagnostic_parameters=['.0001', '.001', '.01', '1e-50'],
                maximum_linear_moment_error='1e-25', original_moment_order_needed_for_loop_C3=4,
                scope='Actual core analytic background coefficient conditional on inherited bounds; unmodulated C3 transfer; conditional linear patch inverse; manufactured diagnostics. Final derivative budget, global extension, waves, infinite sum, smooth force and complete-field blow-up remain open.')


def provenance():
    """Preserve the physical checkpoint and every earlier manuscript audit hash."""
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit',
               'intermediate_decay_audit', 'pulse_moment_audit', 'pulse_stress_audit',
               'post_pulse_stress_audit', 'heat_exterior_audit', 'axis_matching_audit',
               'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds',
               'repair_state_bounds', 'physical_residual_budget')
    checks = {}
    for folder in folders:
        record = json.loads((HERE.parent/folder/'evidence.json').read_text())
        if 'source_sha256' in record:
            inherited = record['source_sha256']
            source_root = HERE.parent/folder
        elif 'source_hashes' in record:
            inherited = record['source_hashes']
            source_root = PROJECT
        else:
            inherited = record.get('provenance', {})
            source_root = None
        if not inherited:
            raise ValueError('Missing inherited hash ledger: '+folder)
        for path, digest in inherited.items():
            root = source_root if source_root is not None else (PROJECT if path.startswith('results/') or path == 'GOAL.md' else HERE.parent)
            checks[folder+':'+path] = hashlib.sha256((root/path).read_bytes()).hexdigest() == digest
    return checks


def hashes():
    names = ('README.md', 'requirements.txt', 'protocol.json', 'transfer.py',
             'system.py', 'majorant.py', 'series.py', 'checks.py', 'audit.py', 'test_audit.py')
    paths = [HERE/name for name in names]
    paths += [PROJECT/'GOAL.md',
              HERE.parent/'physical_residual_budget'/'evidence.json',
              HERE.parent/'axis_core_attachment'/'evidence.json',
              HERE.parent/'angular_c2_transfer'/'evidence.json',
              HERE.parent/'stress_realization_audit'/'evidence.json',
              HERE.parent/'nonlinear_axis_pilot'/'inner.py',
              HERE.parent/'outer_pressure_pilot'/'schedule.py']
    return {str(path.relative_to(PROJECT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def encode(value):
    if hasattr(value, '_mpi_'):
        return dict(lower=mp.nstr(lower(value), 120), upper=mp.nstr(upper(value), 120))
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 80)
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite diagnostic')
        return float(value)
    return value


def protected_growth_ledger():
    return dict(conditional_on_constructed_actual_inner_coefficient=True,
                point='X=1/Lambda, eta=0', physical_path='z=0, r=sqrt(2*q/Lambda), t=1-q',
                strip_loss='rho/16', z='2*C^16*sqrt(3/Lambda)', beta='e*z^2/(rho/16)',
                first_coefficient_bound='M1=(z+beta)*exp(beta)',
                positive_threshold='q0=min(1,(.132/(C*M1))^(1/(2*h)))',
                first_partial_swirl_lower_bound='.132*sqrt(2/Lambda)/C*q^(-A)',
                threshold_materialized=False,
                complete_corrected_field_bound_supplied=False)


def generate():
    mp.mp.dps = 80
    protocol = json.loads((HERE/'protocol.json').read_text())
    if protocol != expected_protocol():
        raise ValueError('Frozen source, domain or scientific scope changed')
    transfer = [angular_transfer(digits) for digits in protocol['interval_digits']]
    core = [inner_operator_certificate(digits) for digits in protocol['interval_digits']]
    constants = norm_ledger()
    reviews = dict(matrices=matrix_review(), sparse_words=sparse_review(),
                   recurrence=recurrence_review(), linear_moments=first_moment_review(),
                   third_implicit=implicit_review())
    old = provenance()
    gates = dict(actual_unmodulated_C3_transfer=all(all(row['checks'].values()) for row in transfer),
                 explicit_actual_core_operator_envelope=all(all(row['checks'].values()) for row in core),
                 third_norm_and_derivative_budget=all(constants['checks'].values()),
                 independent_physical_system_crosscheck=all(reviews['matrices']['checks'].values()),
                 sparse_Picard_majorant=all(reviews['sparse_words']['checks'].values()),
                 manufactured_recurrence_and_axis_slopes=all(reviews['recurrence']['checks'].values()),
                 linear_first_moment_map_and_nonlinear_control=all(reviews['linear_moments']['checks'].values()),
                 raw_third_implicit_derivative=all(reviews['third_implicit']['checks'].values()),
                 inherited_bytes_preserved=bool(old) and all(old.values()),
                 global_and_blowup_scope_preserved=not any(FLAGS[key] for key in
                     ('actual_post_modulation_C3_bound_verified', 'actual_global_first_coefficient_extension_verified',
                      'full_PDE_corrections_verified', 'smooth_force_verified', 'complete_field_energy_verified', 'blowup_verified')))
    return encode(dict(status=STATUS, flags=FLAGS, protocol=protocol, source_hashes=hashes(),
                       unmodulated_C3_transfer=transfer, actual_core_operator_bounds=core,
                       norm_ledger=constants, implicit_derivative_ledger=third_implicit_terms(),
                       analytic_construction=proof_ledger(), protected_partial_growth=protected_growth_ledger(),
                       reviews=reviews, inherited_provenance=old, gates=gates,
                       resolved_obligations=RESOLVED, remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS:
        raise ValueError('An inner analytic coefficient cannot promote the global PDE or blow-up')
    if record.get('resolved_obligations') != RESOLVED or record.get('remaining_obligations') != REMAINING:
        raise ValueError('The unresolved derivative, support or summation obligations changed')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A construction, negative control or provenance check failed')


def compare_records(actual, frozen, location='root'):
    if isinstance(actual, bool) or isinstance(actual, (str, int)) or actual is None:
        if type(actual) is not type(frozen) or actual != frozen:
            raise ValueError('Frozen evidence differs at '+location)
    elif isinstance(actual, float):
        if isinstance(frozen, bool) or not isinstance(frozen, (float, int)) or not math.isclose(actual, frozen, rel_tol=1e-8, abs_tol=2e-10):
            raise ValueError('Frozen numerical evidence differs at '+location)
    elif isinstance(actual, dict):
        if not isinstance(frozen, dict) or set(actual) != set(frozen):
            raise ValueError('Frozen fields differ at '+location)
        for key in actual:
            compare_records(actual[key], frozen[key], location+'.'+key)
    elif isinstance(actual, list):
        if not isinstance(frozen, list) or len(actual) != len(frozen):
            raise ValueError('Frozen samples differ at '+location)
        for i, (a, b) in enumerate(zip(actual, frozen)):
            compare_records(a, b, location+'['+str(i)+']')
    else:
        raise TypeError('Unsupported evidence value')


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
    print(str(len(result['gates']))+' gates pass; conditional actual inner coefficient and linear moment inverse supplied.')
    print('Final derivative bounds, global extension, physical waves, smooth force and complete-field blow-up remain open.')


if __name__ == '__main__':
    main()
