#!/usr/bin/env python3
"""Reproduce local core growth and the frozen-width activation budget."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import mpmath as mp
from diagnostics import algebra_review, budget_review, moment_review, pressure_review, scalar_review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'core-growth-and-quantitative-inner-activation-audited'
FLAGS = dict(unedited_reference_axis_pressure_family_identified=True,
             positive_angular_source_without_division_by_chi_derived=True,
             exact_local_core_shear_and_two_sector_exit_growth_derived=True,
             activation_K_bound_at_same_C_and_width_derived=True,
             complete_local_pressure_and_forward_moments_retained=True,
             local_leading_inner_collar_cone_and_positive_gaps_derived=True,
             completed_global_backward_pressure_identification_verified=False,
             compact_modulation_and_repair_cone_independently_verified=False,
             whole_leading_profile_verified=False,
             actual_giant_schedule_solution_numerically_materialized=False,
             unconditional_uniform_actual_pulse_estimate_verified=False,
             unconditional_complete_residual_improvement_verified=False,
             blowup_verified=False, formal_proof_assistant_certificate=False,
             independent_peer_review_completed=False)
REVIEWS = dict(exact_source_scalar_and_separation_inequalities=algebra_review,
               scalar_comparison_series=scalar_review,
               labelled_pressure_microregion=pressure_review,
               complete_forward_moments_and_C_gauge=moment_review,
               outward_frozen_width_and_growth_budgets=budget_review)
REMAINING = ('audit_compact_modulation_and_exact_repair_cones_with_pressure_memory',
             'identify_completed_global_backward_pressure_using_all_exact_moment_identities',
             'verify_global_coefficient_supports_and_reserved_mean_supports',
             'verify_unconditional_uniform_actual_pulses_across_every_physical_scale',
             'verify_unconditional_complete_residual_gain_including_viscosity_cutoffs_and_interactions')


def protocol():
    return dict(status=STATUS, parent_commit='3108ba7d550f074afc1199188609c898195b525d',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['A.21', 'B.17-B.31', '4.16', '7.1'],
                same_fixed_N_h_Lambda_C_and_reference_datum=True,
                reference_pressure_measure='positive on [0,1], mass<5, ideal atom>=5/2 at theta=1',
                completed_global_pressure_identification_not_implied_by_reference_datum=True,
                j='1e-18', sigma='1e-20', h_range='0<h<=1/100', pressure_range='P>=16',
                core_real_interval='0<=Y<=4.1', core_exit_Y=4,
                angular_source_lower='2.5+.95*Lambda*L*chi',
                two_exit_sectors=['chi>.99', 'chi<=.99'],
                native_exit_growth_lower='2.3', native_collar_growth_lower='9/4',
                activation_K='C^2*exp(.8*B)<C^3',
                unchanged_width_and_kappa0='C^-20', inner_s='1/(10*sqrt(log C))',
                activated_growth_lower='17/8', squared_cone_loss_upper='(32/81)*C^-31*s^6<1/2',
                lambda_squared_lower='exp(-4*Lambda*M)/(680*Lambda*C^2)',
                diagnostic_pressure_atoms=[['1', '2.5'], ['.3', '1.2']],
                diagnostic_gauges=[1, 2, 4], diagnostic_width='.01',
                diagnostic_s='.18', diagnostic_eta='.23', quadrature_order=24,
                diagnostic_precision=140, interval_digits=[100, 160], interval_scales=['P16', 'Md4', 'Md64'],
                manufactured_fixtures_are_actual_core=False,
                fixed_profile_constants_uniform_as_C_tends_to_infinity=False,
                physical_radial_derivatives_uniform_in_width=False,
                local_leading_cone_is_unconditional_global_pulse_or_residual_estimate=False,
                diagnostics_are_not_continuum_or_proof_assistant_certificates=True)


def inherited_hashes():
    parent = HERE.parent/'core_contraction_audit'
    spec = importlib.util.spec_from_file_location('frozen_growth_parent_audit', parent/'audit.py')
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(parent))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    checks = module.inherited_hashes()
    record = json.loads((parent/'evidence.json').read_text())
    for path, expected in record['source_hashes'].items():
        checks['core_contraction_audit:'+path] = hashlib.sha256((PROJECT/path).read_bytes()).hexdigest() == expected
    return checks


def source_hashes():
    paths = [HERE/name for name in ('README.md', '.gitignore', 'requirements.txt', 'protocol.json',
                                   'growth.py', 'diagnostics.py', 'audit.py', 'test_audit.py')]
    paths += [PROJECT/'GOAL.md', HERE.parent/'core_contraction_audit/evidence.json',
              HERE.parent/'inner_collar_activation/evidence.json', HERE.parent/'background_jet_foundation/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != protocol():
        raise ValueError('Frozen growth/activation protocol changed')
    with mp.workdps(140):
        reviews = {name: fn() for name, fn in REVIEWS.items()}
    inherited = inherited_hashes()
    gates = {name: bool(r['checks']) and all(r['checks'].values()) for name, r in reviews.items()}
    gates['inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['local_leading_scope_and_global_obligations_retained'] = not any(FLAGS[k] for k in
        ('completed_global_backward_pressure_identification_verified',
         'compact_modulation_and_repair_cone_independently_verified', 'whole_leading_profile_verified',
         'actual_giant_schedule_solution_numerically_materialized',
         'unconditional_uniform_actual_pulse_estimate_verified',
         'unconditional_complete_residual_improvement_verified', 'blowup_verified'))
    return dict(status=STATUS, protocol=protocol(), flags=dict(FLAGS), reviews=reviews, gates=gates,
                remaining_obligations=list(REMAINING), inherited_provenance=inherited, source_hashes=source_hashes())


def validate(record):
    if record.get('status') != STATUS or record.get('protocol') != protocol() or record.get('flags') != FLAGS:
        raise ValueError('Changed datum, local growth/activation scope or fixed parameters')
    if record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('A global pressure, support, pulse or residual obligation was omitted')
    expected = set(REVIEWS)|{'inherited_bytes_preserved', 'local_leading_scope_and_global_obligations_retained'}
    if set(record.get('gates', {})) != expected or not all(record['gates'].values()):
        raise ValueError('A core growth/activation gate failed')
    if not record.get('inherited_provenance') or not all(record['inherited_provenance'].values()):
        raise ValueError('An inherited source changed')
    if set(record.get('reviews', {})) != set(REVIEWS) or not all(
            r.get('checks') and all(r['checks'].values()) for r in record['reviews'].values()):
        raise ValueError('A complete pressure, source, gauge or log-budget check failed')


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
        if not math.isclose(actual, frozen, rel_tol=3e-9, abs_tol=1e-100):
            raise ValueError('Changed numerical diagnostic at '+path)
    elif isinstance(actual, str) and actual != frozen:
        if path.rsplit('.', 1)[-1] not in ('lower', 'upper') or not isinstance(frozen, str):
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
    print('Local leading inner collar only; compact/global pulses, complete residual gain and blow-up remain unverified.')


if __name__ == '__main__':
    main()
