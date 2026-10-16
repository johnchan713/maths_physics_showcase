#!/usr/bin/env python3
"""Reproduce the conditional inner activation and exact-seed obstruction."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from review import (activation_review, endpoint_review, obstruction_review,
                    primitive_review, quotient_review, refinement_review)

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'inner-flat-activation-and-conditional-cone-audited'
FLAGS = dict(width_rescaled_flat_Volterra_identity_derived=True,
             full_forward_moments_and_pressure_difference_retained=True,
             conditional_cubic_direction_and_sixth_order_cone_loss_derived=True,
             common_parameter_collar_propagation_relative_to_exact_seed=True,
             nonzero_core_defect_obstruction_explicit=True,
             actual_core_contraction_independently_reverified=False,
             actual_seed_growth_constants_and_K_numerically_certified=False,
             complete_inner_and_compact_foundation_verified=False,
             global_exact_moment_identification_independently_reverified=False,
             unconditional_uniform_actual_pulse_estimate_verified=False,
             unconditional_complete_residual_improvement_verified=False,
             blowup_verified=False, formal_proof_assistant_certificate=False)
REVIEWS = dict(flat_quotient_and_endpoint_jets=quotient_review,
               five_forward_primitive_derivatives=primitive_review,
               complete_activation_difference_and_eta_jets=activation_review,
               endpoint_coefficients_and_explicit_collar=endpoint_review,
               tighter_positive_quadrature=refinement_review,
               baseline_defect_and_failed_shortcuts=obstruction_review)
REMAINING = ('independently_verify_actual_nonlinear_core_contraction_and_exact_stress_zero_identity',
             'certify_actual_seed_positive_gaps_and_normalized_error_bound_K_at_same_N_h',
             'independently_audit_compact_modulation_repair_and_global_exact_moments',
             'verify_reserved_mean_supports_and_actual_global_assembly',
             'verify_unconditional_uniform_actual_pulses_and_complete_residual_gain')


def protocol():
    return dict(status=STATUS, parent_commit='da0b3fe4a297f0367b75d095fac029103aa0153c',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['4.7-4.11', '4.15-4.16', '4.20-4.23', 'B.26-B.31', '7.1'],
                same_fixed_N_and_h=True, theorem_seed_is_exact_zero_stress_core=True,
                theorem_dimensionless_domain='0<=s=y/t1<=1/4, 0<t1<=tbar, 0<=kappa0<=1/2',
                diagnostic_width_max='.01', diagnostic_h='.005', diagnostic_precision=80,
                quadrature_orders=[24, 40], primitive_seed_Xa=2,
                manufactured_seed_exact_axis_solution=False,
                diagnostic_cone_is_activated_excess_not_full_seed_stress=True,
                field_factor='t1*ea*s^3', moment_factor='t1^2*ea*s^6',
                normalized_direction_error='O(t1*s^3)', squared_cone_loss='O(t1^2*s^6)',
                physical_radial_derivatives_uniform_in_t1=False,
                actual_core_tolerance_cannot_replace_exact_identity=True,
                diagnostics_are_not_continuum_proofs=True)


def inherited_hashes():
    # Re-evaluate every earlier source ledger, not merely its stored booleans.
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit', 'intermediate_decay_audit',
               'pulse_moment_audit', 'pulse_stress_audit', 'post_pulse_stress_audit', 'heat_exterior_audit',
               'axis_matching_audit', 'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds', 'repair_state_bounds',
               'physical_residual_budget', 'first_correction_construction', 'higher_jet_budget',
               'annular_first_correction', 'background_summation', 'physical_wave_interface',
               'actual_band_pulse_bounds', 'compact_mean_corrections', 'residual_cycle_closure',
               'full_potential_realization', 'localized_force_extension', 'global_solution_comparison',
               'auxiliary_geometry_audit', 'geometric_pulse_replay', 'background_jet_foundation', 'outer_edge_stability')
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
            raise ValueError('Missing inherited source ledger: '+folder)
        for path, expected in hashes.items():
            routed = root if root is not None else (PROJECT if path.startswith('results/') or path == 'GOAL.md' else HERE.parent)
            checks[folder+':'+path] = hashlib.sha256((routed/path).read_bytes()).hexdigest() == expected
    return checks


def source_hashes():
    names = ('README.md', '.gitignore', 'requirements.txt', 'protocol.json', 'collar.py',
             'review.py', 'audit.py', 'test_audit.py')
    paths = [HERE/name for name in names]+[PROJECT/'GOAL.md',
             HERE.parent/'outer_edge_stability/evidence.json',
             HERE.parent/'axis_core_attachment/evidence.json',
             HERE.parent/'background_jet_foundation/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != protocol():
        raise ValueError('Frozen activation protocol changed')
    with mp.workdps(80):
        reviews = {name: fn() for name, fn in REVIEWS.items()}
    inherited = inherited_hashes()
    gates = {name: bool(r['checks']) and all(r['checks'].values()) for name, r in reviews.items()}
    gates['inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['actual_seed_and_global_scope_not_promoted'] = not any(FLAGS[k] for k in
        ('actual_core_contraction_independently_reverified', 'actual_seed_growth_constants_and_K_numerically_certified',
         'complete_inner_and_compact_foundation_verified', 'global_exact_moment_identification_independently_reverified',
         'unconditional_uniform_actual_pulse_estimate_verified',
         'unconditional_complete_residual_improvement_verified', 'blowup_verified'))
    return dict(status=STATUS, protocol=protocol(), flags=dict(FLAGS), reviews=reviews, gates=gates,
                remaining_obligations=list(REMAINING), inherited_provenance=inherited, source_hashes=source_hashes())


def validate(record):
    if record.get('status') != STATUS or record.get('protocol') != protocol() or record.get('flags') != FLAGS:
        raise ValueError('Changed activation scope, normalization or exact-seed premise')
    if record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('An actual-core or global obligation was omitted')
    expected = set(REVIEWS)|{'inherited_bytes_preserved', 'actual_seed_and_global_scope_not_promoted'}
    if set(record.get('gates', {})) != expected or not all(record['gates'].values()):
        raise ValueError('A conditional inner activation gate failed')
    if not record.get('inherited_provenance') or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited bytes changed')
    if set(record.get('reviews', {})) != set(REVIEWS) or not all(
            r.get('checks') and all(r['checks'].values()) for r in record['reviews'].values()):
        raise ValueError('A full-memory or failed-shortcut check is absent')


def compare_records(actual, frozen, path='root'):
    if isinstance(actual, dict):
        if not isinstance(frozen, dict) or set(actual) != set(frozen):
            raise ValueError('Changed record keys at '+path)
        for key in actual:
            compare_records(actual[key], frozen[key], path+'.'+key)
    elif isinstance(actual, list):
        if not isinstance(frozen, list) or len(actual) != len(frozen):
            raise ValueError('Changed sequence at '+path)
        for i, value in enumerate(actual):
            compare_records(value, frozen[i], path+'.'+str(i))
    elif isinstance(actual, float):
        if type(frozen) is not float or not math.isfinite(actual) or not math.isfinite(frozen):
            raise ValueError('Nonfinite or mistyped result at '+path)
        # Tiny errors may vary by precision/platform. Do not allow a nonzero
        # coefficient or defect to be erased with an absolute tolerance.
        is_error = any(k in path for k in ('relative_error', 'derivative_error', 'jet_scaled_error',
                                          'identity_error', 'difference_error'))
        if not math.isclose(actual, frozen, rel_tol=3e-9, abs_tol=1e-12 if is_error else 1e-300):
            raise ValueError('Changed numerical result at '+path)
    elif type(actual) is not type(frozen) or actual != frozen:
        raise ValueError('Changed exact result at '+path)


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
    args.output.write_text(json.dumps(record, sort_keys=True, indent=2)+'\n')
    print(STATUS)
    print(f"{len(record['gates'])} gates pass; {len(record['inherited_provenance'])} inherited hashes agree.")
    print('Exact actual core, unconditional pulses, complete residual gain and blow-up remain unverified.')


if __name__ == '__main__':
    main()
