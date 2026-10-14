#!/usr/bin/env python3
"""Reproduce scoped higher jets without promoting singular residuals to blow-up."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import mpmath as mp

from higher_bounds import (HERE, PROJECT, INPUT_POWER, lower, upper, cutoff_certificate,
                           transfer_C4, compact_exponents, scale_checks, loop_ledger,
                           final_ledger, fourth_implicit_terms)
from higher_checks import (variance_review, root_review, phase_review, modulation_review,
                           fourth_implicit_review, repair_C3_review, moment_review,
                           derivative_gap_review, reference_review)
import sys
sys.path.insert(0, str(HERE))

STATUS = 'conditional-finalized-compact-background-jets-bounded'
FLAGS = dict(conditional_actual_original_C4_transfer_supplied=True,
             conditional_actual_mixed_input_envelope_supplied=True,
             conditional_actual_loop_J3_bound_supplied=True,
             conditional_actual_final_repair_C3_bound_supplied=True,
             conditional_actual_compact_final_K_supplied=True,
             accepted_frequency_and_actual_profile_unchanged=True,
             relies_on_inherited_analytic_profile_and_value_bounds=True,
             actual_high_frequency_profile_numerically_resolved=False,
             actual_global_first_coefficient_extension_verified=False,
             full_physical_stress_wave_construction_verified=False,
             full_physical_residual_budget_verified=False,
             physical_background_residual_decay_verified=False,
             full_PDE_corrections_verified=False,
             smooth_force_verified=False, complete_field_energy_verified=False,
             blowup_verified=False, independent_foundation_review_completed=False,
             formal_proof_assistant_certificate=False)
RESOLVED = ['conditional_original_C4_and_loop_J3_bounds_on_selected_compact_annulus',
            'conditional_final_repair_C3_bound_at_same_frequency',
            'conditional_ordinary_C2_annular_profile_constant_K_H64']
REMAINING = ['actual_first_coefficient_annular_extension_and_five_moment_support',
             'full_wave_covariance_and_all_nonlinear_interactions',
             'complete_physical_residual_budget_and_localization',
             'all_order_PDE_correction_summation',
             'smooth_compact_forcing_through_T',
             'complete_field_blowup_lower_bound_and_energy',
             'independent_review_of_inherited_continuum_estimates']


def expected_protocol():
    return dict(status=STATUS, parent_commit='007e62724850fbe98707d5c50578e47b35768c4d',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['4.9', '4.16', 'B.22', 'B.25-B.26', 'C.5', 'C.9-C.16'],
                interval_digits=[80, 110], fixed_N='1+floor(H^32)', original_A='C^4096',
                higher_derivative_envelope='A^256', loop_joint_order=3,
                original_mixed_field_and_moment_order=4, final_angular_moment_order=3,
                ordinary_profile_jet_order=2, actual_annular_K='H^64', fixed_viscosity=1,
                physical_q_range='0<q<=1', domain='X in [Xa*exp(tc/2), rhop*exp(5)], eta in [-1,1]',
                step_derivative_maximum_order=5, step_derivative_upper_bound='1e100',
                maximum_multiprecision_identity_error='1e-55',
                maximum_legacy_binary64_moment_error='2e-14',
                maximum_nonlinear_scaled_identity_error='1e-50',
                manufactured_diagnostic_frequencies=[8, 32, 128],
                counterexample_frequencies=[8, 32, 128, 512],
                actual_high_frequency_profile_numerically_resolved=False,
                scope='Conditional actual compact-annular input C4, loop J3, final-repair C3 and ordinary background C2 bounds at the unchanged N. Singular background residual bounds are explicit; wave stress, global first correction, all-order summation, smooth forcing, energy and blow-up remain open.')


def provenance():
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit',
               'intermediate_decay_audit', 'pulse_moment_audit', 'pulse_stress_audit',
               'post_pulse_stress_audit', 'heat_exterior_audit', 'axis_matching_audit',
               'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds',
               'repair_state_bounds', 'physical_residual_budget', 'first_correction_construction')
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
    paths = [HERE/name for name in ('README.md', 'requirements.txt', 'protocol.json',
                                    'higher_bounds.py', 'higher_jets.py', 'higher_checks.py',
                                    'audit.py', 'test_audit.py')]
    paths += [PROJECT/'GOAL.md']
    for folder in ('first_correction_construction', 'axis_core_attachment', 'compact_jet_envelope',
                   'loop_modulation_bounds', 'repair_state_bounds', 'physical_residual_budget'):
        paths.append(HERE.parent/folder/'evidence.json')
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
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    return value


def generate():
    mp.mp.dps = 80
    protocol = json.loads((HERE/'protocol.json').read_text())
    if protocol != expected_protocol():
        raise ValueError('Frozen domain, source or derivative scope changed')
    cutoff = [cutoff_certificate(d) for d in protocol['interval_digits']]
    transfer = [transfer_C4(d) for d in protocol['interval_digits']]
    scales = [scale_checks(d) for d in protocol['interval_digits']]
    inputs, loops, final = compact_exponents(), loop_ledger(), final_ledger()
    reviews = dict(reference=reference_review(), variance=variance_review(), implicit_root=root_review(), moving_phase=phase_review(),
                   modulation=modulation_review(), fourth_implicit=fourth_implicit_review(),
                   final_repair=repair_C3_review(), moments=moment_review(), gap=derivative_gap_review())
    inherited = provenance()
    gates = dict(uniform_cutoff_through_order_five=all(all(x['checks'].values()) for x in cutoff),
                 derived_actual_C4_reference_and_transfer=all(all(x['checks'].values()) for x in transfer),
                 independent_reference_source_integral=all(reviews['reference']['checks'].values()),
                 mixed_input_and_denominator_costs=all(inputs['checks'].values()),
                 unchanged_scale_comparisons=all(all(x['checks'].values()) for x in scales),
                 regular_variance_and_signed_root=all(reviews['variance']['checks'].values()),
                 third_implicit_root_and_boundary=all(reviews['implicit_root']['checks'].values()),
                 third_moving_phase_and_negative_controls=all(reviews['moving_phase']['checks'].values()),
                 fixed_N_and_exponential_chain_controls=all(reviews['modulation']['checks'].values()),
                 fourth_joining_derivative_terms=all(reviews['fourth_implicit']['checks'].values()),
                 final_C3_repair_and_normalization=all(reviews['final_repair']['checks'].values()),
                 five_original_C4_moments=all(reviews['moments']['checks'].values()),
                 derivative_gap_controls=all(reviews['gap']['checks'].values()),
                 analytic_loop_and_final_budget_ledger=all(loops['checks'].values()) and all(final['checks'].values()),
                 inherited_bytes_preserved=bool(inherited) and all(inherited.values()),
                 singular_background_bounds_do_not_promote_blowup=not any(FLAGS[key] for key in
                    ('full_physical_residual_budget_verified', 'physical_background_residual_decay_verified',
                     'actual_global_first_coefficient_extension_verified', 'full_PDE_corrections_verified',
                     'smooth_force_verified', 'complete_field_energy_verified', 'blowup_verified')))
    return encode(dict(status=STATUS, protocol=protocol, flags=FLAGS, cutoff=cutoff, C4_transfer=transfer,
                       input_exponents=inputs, scale_comparisons=scales, loop_derivatives=loops,
                       final_budget=final, fourth_implicit= fourth_implicit_terms(), reviews=reviews,
                       inherited_provenance=inherited, source_hashes=hashes(), gates=gates,
                       resolved_obligations=RESOLVED, remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS:
        raise ValueError('Compact jets and singular upper bounds cannot promote the global PDE or blow-up')
    if record.get('resolved_obligations') != RESOLVED or record.get('remaining_obligations') != REMAINING:
        raise ValueError('Unclosed extension, wave, summation or forcing obligations changed')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A higher-jet estimate, diagnostic or provenance check failed')


def compare_records(actual, frozen, location='root'):
    if isinstance(actual, dict):
        if not isinstance(frozen, dict) or set(actual) != set(frozen):
            raise ValueError('Frozen fields differ at '+location)
        for key in actual:
            compare_records(actual[key], frozen[key], location+'.'+key)
    elif isinstance(actual, list):
        if not isinstance(frozen, list) or len(actual) != len(frozen):
            raise ValueError('Frozen samples differ at '+location)
        for i, (a, b) in enumerate(zip(actual, frozen)):
            compare_records(a, b, location+'['+str(i)+']')
    elif isinstance(actual, str) and isinstance(frozen, str) and location.startswith('root.reviews.'):
        try:
            a, b = mp.mpf(actual), mp.mpf(frozen)
        except ValueError:
            if actual != frozen:
                raise ValueError('Frozen review label differs at '+location)
        else:
            if not mp.isfinite(a) or not mp.isfinite(b) or abs(a-b) > max(mp.mpf('1e-65'), mp.mpf('1e-12')*max(abs(a), abs(b))):
                raise ValueError('Frozen numerical review differs at '+location)
    elif type(actual) is not type(frozen) or actual != frozen:
        raise ValueError('Frozen evidence differs at '+location)


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
    print(str(len(result['gates']))+' gates pass; conditional compact finalized K=H^64 at the unchanged N.')
    print('Residual bounds are singular; global coefficient extension, physical waves, force and blow-up remain open.')


if __name__ == '__main__':
    main()
