#!/usr/bin/env python3
"""Reproduce the scoped background foundation audit and preserve earlier evidence."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from diagnostics import (axis_review, chart_review, derivative_review,
                         interval_matrix_review, matrix_quadrature_review,
                         picard_review, scale_review)

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'background-jet-construction-audited-relative-to-leading-profile'
FLAGS = dict(scale_derivative_conversion_and_diagonal_constants_explicit=True,
             regular_axis_quotients_and_cutoff_curl_retained=True,
             all_order_common_collar_argument_supplied=True,
             actual_bump_matrix_continuum_inverse_below_1000=True,
             normalized_background_jet_bound_derived_relative_to_leading_certificate=True,
             common_band_domain_includes_one_sided_boundary=True,
             actual_leading_profile_or_cutoff_scales_numerically_evaluated=False,
             full_leading_profile_foundation_independently_verified=False,
             strict_actual_pulse_gap_and_weighted_edge_certificate_verified=False,
             unconditional_uniform_pulse_estimate_verified=False,
             unconditional_complete_residual_improvement_verified=False,
             blowup_verified=False, formal_proof_assistant_certificate=False)
REVIEWS = dict(scale_derivatives=derivative_review, regular_axis=axis_review,
               continuum_moment_inverse=interval_matrix_review,
               actual_bump_quadrature=matrix_quadrature_review,
               whole_interval_cutoff_schedule=scale_review,
               common_collar_picard=picard_review, band_pullback=chart_review)
REMAINING = ('verify_finalized_leading_profile_regular_analytic_collar_and_patch_exterior_identities_at_same_N_h',
             'verify_actual_strict_growth_cone_gaps_and_weighted_T0_T1_edges',
             'audit_reserved_physical_mean_pressure_support_and_selected_localization_potentials',
             'independent_complete_residual_and_global_assembly_review_before_unconditional_claims')


def protocol():
    return dict(status=STATUS, parent_commit='94e1f3376d6813376d5369bba9c3b949a9802dd4',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                actual_N_h='same fixed inherited N and h; neither reselected',
                derivative_convention='Euler q*d_q at fixed X,eta; Stirling conversion from q^j*d_q^j',
                input='leading-profile regularity, common analytic collar, reserved patch and exterior identities',
                precision=80, interval_precision=60,
                exact_derivative_orders=list(range(11)),
                diagnostic_patch_exponents=['0.005', '1e-30'],
                actual_patch_exponent_range='0<ell_patch<=1/100',
                bump_width='1/10', U_bump_centers=['1', '2'], E_bump_centers=['1/2', '3/2', '5/2'],
                common_collar='0<=X<=Xa*exp(tc/4); keep through Xa*exp(tc/8)',
                normalized_diagonal='C_(n,m)*q^(n*h)<=2^-n for every m<=n and all q<=a_n^-1',
                radial_cutoff_jet_order='m+1', sampled_diagnostics_are_not_continuum_proofs=True)


def inherited_hashes():
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit', 'intermediate_decay_audit',
               'pulse_moment_audit', 'pulse_stress_audit', 'post_pulse_stress_audit', 'heat_exterior_audit',
               'axis_matching_audit', 'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds', 'repair_state_bounds',
               'physical_residual_budget', 'first_correction_construction', 'higher_jet_budget',
               'annular_first_correction', 'background_summation', 'physical_wave_interface',
               'actual_band_pulse_bounds', 'compact_mean_corrections', 'residual_cycle_closure',
               'full_potential_realization', 'localized_force_extension', 'global_solution_comparison',
               'auxiliary_geometry_audit', 'geometric_pulse_replay')
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
        for path, expected in hashes.items():
            routed = root if root is not None else (PROJECT if path.startswith('results/') or path == 'GOAL.md' else HERE.parent)
            checks[folder+':'+path] = hashlib.sha256((routed/path).read_bytes()).hexdigest() == expected
    return checks


def source_hashes():
    names = ('README.md', '.gitignore', 'requirements.txt', 'protocol.json', 'foundation.py',
             'diagnostics.py', 'audit.py', 'test_audit.py')
    paths = [HERE/name for name in names]+[PROJECT/'GOAL.md', PROJECT/'EXTERNAL_REVIEW_2026-10-15.md',
                                         HERE.parent/'background_summation/evidence.json',
                                         HERE.parent/'geometric_pulse_replay/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != protocol():
        raise ValueError('Protocol changed')
    with mp.workdps(80):
        reviews = {name: fn() for name, fn in REVIEWS.items()}
    inherited = inherited_hashes()
    gates = {name: all(review['checks'].values()) for name, review in reviews.items()}
    gates['inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['leading_profile_and_pulse_gap_still_open'] = not any(FLAGS[name] for name in
        ('full_leading_profile_foundation_independently_verified',
         'strict_actual_pulse_gap_and_weighted_edge_certificate_verified',
         'unconditional_uniform_pulse_estimate_verified',
         'unconditional_complete_residual_improvement_verified', 'blowup_verified'))
    return dict(status=STATUS, protocol=protocol(), flags=FLAGS, reviews=reviews,
                gates=gates, remaining_obligations=list(REMAINING),
                inherited_provenance=inherited, source_hashes=source_hashes())


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS or record.get('protocol') != protocol():
        raise ValueError('Scope, derivative convention or construction inputs changed')
    if record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('Leading certificate, pulse gaps or global review omitted')
    if set(record.get('gates', {})) != set(REVIEWS)|{'inherited_bytes_preserved', 'leading_profile_and_pulse_gap_still_open'}:
        raise ValueError('Missing or additional audit gate')
    if not all(record['gates'].values()):
        raise ValueError('A required audit gate failed')
    if not record.get('inherited_provenance') or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited source changed')
    if set(record.get('reviews', {})) != set(REVIEWS):
        raise ValueError('Missing review')
    if not all(all(r.get('checks', {}).values()) and r.get('checks') for r in record['reviews'].values()):
        raise ValueError('Failed diagnostic')


def compare_records(actual, frozen, path='root'):
    if isinstance(actual, dict):
        if not isinstance(frozen, dict) or set(actual) != set(frozen):
            raise ValueError('Changed keys at '+path)
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
        tolerance = 1e-45 if 'derivative_error' in path else (1e-15 if 'refinement_error' in path else 1e-12)
        if not math.isclose(actual, frozen, rel_tol=1e-11, abs_tol=tolerance):
            raise ValueError('Changed numerical diagnostic at '+path)
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
    print('Leading profile, strict pulse gaps, unconditional residual and blow-up remain unverified.')


if __name__ == '__main__':
    main()
