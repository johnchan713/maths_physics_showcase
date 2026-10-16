#!/usr/bin/env python3
"""Reproduce the core contraction audit with its exact pressure-family premise."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import mpmath as mp
from checks import (budget_review, failure_review, formal_review, geometry_review,
                    inverse_review, source_review, weight_review)

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'core-contraction-and-exact-stress-audited-for-pressure-family'
FLAGS = dict(all_index_integrated_derivative_and_product_bounds_derived=True,
             nonlinear_contraction_for_stated_exact_pressure_family_audited=True,
             universal_contraction_ceiling_below_1e_minus50_derived=True,
             exact_core_pressure_balance_and_zero_leading_stress_derived=True,
             common_analytic_axis_neighborhood_derived=True,
             actual_global_pressure_family_identification_independently_verified=False,
             actual_core_positive_shear_growth_and_activation_K_verified=False,
             whole_leading_profile_verified=False,
             actual_giant_schedule_solution_numerically_materialized=False,
             unconditional_uniform_actual_pulse_estimate_verified=False,
             unconditional_complete_residual_improvement_verified=False,
             blowup_verified=False, formal_proof_assistant_certificate=False,
             independent_peer_review_completed=False)
REVIEWS = dict(coefficient_weights_and_derivative_allocation=weight_review,
               pressure_family_complex_neighborhood=geometry_review,
               complete_source_and_evaluation_budgets=source_review,
               order_one_angular_inverse=inverse_review,
               original_core_equations_and_forward_stress_jets=formal_review,
               independent_parameter_arithmetic=budget_review,
               failed_real_only_and_full_residual_shortcuts=failure_review)
REMAINING = ('independently_identify_actual_exterior_pressure_with_the_exact_analytic_family',
             'verify_actual_positive_shear_endpoint_growth_and_activation_K_at_same_N_h',
             'audit_compact_modulation_and_repair_with_complete_pressure_memory',
             'verify_global_exact_moments_coefficient_support_and_reserved_mean_support',
             'verify_unconditional_actual_pulses_and_complete_residual_improvement')


def protocol():
    return dict(status=STATUS, parent_commit='5ac315d13e54797b3a8c692a3d465e41d55642ea',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['4.9-4.13', 'B.1-B.16'], same_fixed_N_and_h=True,
                theorem_pressure_family='exact positive measure on [0,1], mass<5; P>=16',
                exact_pressure_identification_is_an_input=True,
                j='1e-18', sigma='1e-20', rho='sigma/2048', h_range='0<h<=1/100',
                coefficient_space='all radial degrees and all eta derivatives',
                universal_contraction_ceiling='1e-50', normalized_jet_orders=[2, 6],
                common_eta_radius='3*rho/8 on |Y|<=5', core_real_interval='0<=Y<=4.1',
                same_selected_Lambda_and_C_formulas=True,
                interval_scales=['P16', 'Md4', 'Md64'], interval_digits=[80, 120],
                formal_diagnostic_degrees=[2, 3, 4], formal_diagnostic_eta_order=2,
                formal_diagnostics_are_actual_core=False,
                zero_leading_stress_is_zero_complete_residual=False,
                finite_diagnostics_are_not_all_index_or_existence_proofs=True)


def inherited_hashes():
    # Reuse the frozen routing function, then verify its own full source ledger.
    parent = HERE.parent/'inner_collar_activation'
    spec = importlib.util.spec_from_file_location('frozen_inner_activation_audit', parent/'audit.py')
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(parent))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    checks = module.inherited_hashes()
    record = json.loads((parent/'evidence.json').read_text())
    for path, expected in record['source_hashes'].items():
        checks['inner_collar_activation:'+path] = hashlib.sha256((PROJECT/path).read_bytes()).hexdigest() == expected
    return checks


def source_hashes():
    paths = [HERE/name for name in ('README.md', '.gitignore', 'requirements.txt', 'protocol.json',
                                   'core.py', 'checks.py', 'audit.py', 'test_audit.py')]
    paths += [PROJECT/'GOAL.md', HERE.parent/'inner_collar_activation/evidence.json',
              HERE.parent/'axis_core_attachment/evidence.json', HERE.parent/'background_jet_foundation/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != protocol():
        raise ValueError('Frozen core-contraction protocol changed')
    with mp.workdps(140):
        reviews = {name: fn() for name, fn in REVIEWS.items()}
    inherited = inherited_hashes()
    gates = {name: bool(r['checks']) and all(r['checks'].values()) for name, r in reviews.items()}
    gates['inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['pressure_premise_growth_and_global_limits_retained'] = not any(FLAGS[k] for k in
        ('actual_global_pressure_family_identification_independently_verified',
         'actual_core_positive_shear_growth_and_activation_K_verified', 'whole_leading_profile_verified',
         'actual_giant_schedule_solution_numerically_materialized',
         'unconditional_uniform_actual_pulse_estimate_verified',
         'unconditional_complete_residual_improvement_verified', 'blowup_verified'))
    return dict(status=STATUS, protocol=protocol(), flags=dict(FLAGS), reviews=reviews, gates=gates,
                remaining_obligations=list(REMAINING), inherited_provenance=inherited, source_hashes=source_hashes())


def validate(record):
    if record.get('status') != STATUS or record.get('protocol') != protocol() or record.get('flags') != FLAGS:
        raise ValueError('Core theorem scope, exact pressure premise or normalization changed')
    if record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('An actual growth, pressure or global obligation was omitted')
    expected = set(REVIEWS)|{'inherited_bytes_preserved', 'pressure_premise_growth_and_global_limits_retained'}
    if set(record.get('gates', {})) != expected or not all(record['gates'].values()):
        raise ValueError('A core-contraction audit gate failed')
    if not record.get('inherited_provenance') or not all(record['inherited_provenance'].values()):
        raise ValueError('An inherited source changed')
    if set(record.get('reviews', {})) != set(REVIEWS) or not all(
            r.get('checks') and all(r['checks'].values()) for r in record['reviews'].values()):
        raise ValueError('A source, inverse, exact-identity or failure-control check failed')


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
            raise ValueError('Nonfinite or mistyped diagnostic at '+path)
        if not math.isclose(actual, frozen, rel_tol=3e-9, abs_tol=1e-60):
            raise ValueError('Changed numerical diagnostic at '+path)
    elif isinstance(actual, str) and actual != frozen:
        numeric = path.rsplit('.', 1)[-1] in ('lower', 'upper', 'log_contraction_upper',
                    'log_normalized_jet_error_upper', 'pressure_over_P2_upper', 'Z_bound_ratio_upper')
        if not numeric or not isinstance(frozen, str):
            raise ValueError('Changed exact string at '+path)
        with mp.workdps(160):
            a, b = mp.mpf(actual), mp.mpf(frozen)
            if not mp.isfinite(a) or not mp.isfinite(b) or abs(a-b) > mp.mpf('1e-55')*max(1, abs(a), abs(b)):
                raise ValueError('Changed interval diagnostic at '+path)
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
    print('Actual pressure identification, growth, unconditional pulses and complete residual gain remain unverified.')


if __name__ == '__main__':
    main()
