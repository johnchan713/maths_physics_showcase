#!/usr/bin/env python3
"""Conditional full-potential sum and flat local residual, not a blow-up proof."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from realization_bounds import (all_stage_order_review, comparison_budget,
                                core_energy_review, scale_choice_review)
from cutoff_checks import (cutoff_geometry_review, curl_review,
                           divergent_series_review, residual_review)

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'conditional-full-potential-realization-and-flat-local-residual'
FLAGS = dict(conditional_full_wave_and_mean_potential_sum_supplied=True,
             every_curl_cutoff_commutator_retained=True,
             conditional_fixed_stage_comparison_and_all_jet_residual_flatness_supplied=True,
             conditional_complete_local_field_core_lower_bound_preserved=True,
             conditional_local_kinetic_energy_bound_supplied=True,
             actual_global_profile_and_field_numerically_evaluated=False,
             actual_infinite_stage_coefficients_or_cutoff_scales_evaluated=False,
             inherited_foundation_and_global_auxiliary_geometry_independently_verified=False,
             closed_time_auxiliary_atlas_and_one_sided_potential_jets_independently_verified=False,
             raw_uncut_infinite_series_convergence_verified=False,
             unconditional_complete_residual_construction_verified=False,
             admissible_smooth_compact_force_verified=False,
             whole_space_field_energy_verified=False,
             global_solution_exclusion_argument_verified=False,
             blowup_verified=False, diagnostic_success_is_not_a_continuum_proof=True,
             formal_proof_assistant_certificate=False)
REVIEW_NAMES = ('all_stage_potential_and_fixed_jet_orders', 'actual_q_cutoff_geometry',
                'complete_curl_cutoff_product_rule', 'complete_cutoff_residual',
                'whole_interval_scale_choice', 'fixed_stage_flat_residual_comparison',
                'protected_core_and_local_energy', 'raw_series_and_quantifier_controls')
GATE_NAMES = REVIEW_NAMES+('all_inherited_bytes_preserved', 'conditional_local_result_not_promoted')
REMAINING = ('closed_time_auxiliary_atlas_and_one_sided_potential_jets_away_from_singular_point',
             'complete_physical_localization_initialization_and_smooth_compact_force_extension',
             'whole_space_energy_and_global_solution_exclusion_argument',
             'independent_inherited_profile_background_and_global_auxiliary_geometry_review')


def expected_protocol():
    return dict(status=STATUS, parent_commit='d80f6ff1b45b6111a400f04db463a9e3afd84ccf',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['5.28-5.33', '9.15-9.21', 'Lemma 9.8', 'Proposition 9.9'],
                actual_profile_parameters='unchanged inherited N, h, positive lambda and fixed maps',
                increment_stage_indexing='j-1 to j; B_(j-1)=3/5+j/10; j>=1',
                fixed_kappa='1/100000', increment_gain='g_j=h*j/10',
                tail_gain='g_j/2=h*j/20', fixed_cutoff_loss='L_m=2*A+2*m+2',
                finite_total_field_loss='H_m=L_m+1',
                finite_residual_loss='K_m=2*A+3/2+2*m',
                diagnostic_frequencies=[11, 19, 29], diagnostic_cutoff_arguments=[.73, .81],
                diagnostic_jet_points=[[1.12, -.4, -.19, .36], [1.3, .27, .2, .47], [1.66, -.9, .31, .64]],
                diagnostic_h='1/200', diagnostic_q_values=[.2, .02, .002],
                diagnostic_float_relative_tolerance=2e-10,
                diagnostic_data='manufactured full wave/mean potentials; actual implicit q; synthetic stage constants',
                scope='Conditional full-potential realization, exact divergence, all-jet flat local residual, preserved core growth and bounded local kinetic energy. Endpoint atlas review, admissible forcing, whole-space conclusions and independent foundation review remain open.')


def proof_ledger():
    return dict(representatives='Delta u_j=curl(Aw_j+Psi_j e_theta)+b_j e_theta; partial_theta b_j=0',
                pressure='complete final pressure differences telescope across every cycle',
                cutoffs='chi(a_j*q), with a_(j+1)>=2*a_j and derivatives supported on 1/2<=a_j*q<=1',
                curl_product_rule='curl(chi*A)=chi*curl(A)+grad(chi) cross A',
                mean_radial_commutator='-a_j*chi_prime(a_j*q)*q_z*Psi_j',
                cutoff_derivative_bound='|partial_z^a partial_t^b chi(a_j*q)|<=C_ab*q^(-D*a-b), independent of a_j',
                stage_constants_may_grow_arbitrarily_fast=True,
                raw_series_convergence_not_used=True,
                tail='|U_loc-U^[J]|_m<=2^(-J)*q^(h*(J+1)/20-L_m), q<=1/(2*a_J), J>=max(1,m)',
                comparison='one fixed J chosen for m,N BEFORE q tends to zero',
                residual_identity='e_t-Delta e+grad pi+u_J.grad e+e.grad u_J+e.grad e',
                sufficient_tail_condition='h*(J+1)/20>=N+1+L_(m+2)+H_(m+2)',
                sufficient_finite_condition='h*(7/10+J/10)-K_m>=N+1',
                finite_flat_errors_compared_not_independently_summed=True,
                real_smooth_divergence_free_for_q_positive=True,
                unchanged_inner_core_and_heat_exterior=True,
                flatness_in_cartesian_distance='q<=2*tau+(2*z^2)^(1/(1-2*h))<=C*(|x|+|t-1|)',
                full_one_sided_endpoint_requires_separate_finite_atlas_review=True,
                normalized_annular_perturbation_bound='|u_loc-u_base|<=C*q^(-A+h/4)',
                conditional_inner_swirl='tau^(-A)*(e0+O(tau^(2*h))), e0>0 inherited',
                local_energy='radial energy<=C*q^(-2*h); axial integral |z|^(-2*h/D) is finite since D-2*h>0',
                local_energy_not_a_whole_space_admissible_solution=True)


def inherited_provenance():
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit',
               'intermediate_decay_audit', 'pulse_moment_audit', 'pulse_stress_audit',
               'post_pulse_stress_audit', 'heat_exterior_audit', 'axis_matching_audit',
               'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds',
               'repair_state_bounds', 'physical_residual_budget', 'first_correction_construction',
               'higher_jet_budget', 'annular_first_correction', 'background_summation',
               'physical_wave_interface', 'actual_band_pulse_bounds', 'compact_mean_corrections',
               'residual_cycle_closure')
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
        for path, digest in hashes.items():
            routed = root if root is not None else (PROJECT if path.startswith('results/') or path == 'GOAL.md' else HERE.parent)
            checks[folder+':'+path] = hashlib.sha256((routed/path).read_bytes()).hexdigest() == digest
    return checks


def source_hashes():
    names = ('README.md', '.gitignore', 'requirements.txt', 'protocol.json',
             'realization_bounds.py', 'cutoff_checks.py', 'audit.py', 'test_audit.py')
    paths = [HERE/name for name in names]+[PROJECT/'GOAL.md',
        HERE.parent/'residual_cycle_closure/evidence.json',
        HERE.parent/'background_summation/evidence.json',
        HERE.parent/'physical_wave_interface/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def encode(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite diagnostic')
        return float(value)
    return value


def comparison_review():
    rows = [comparison_budget(m, N, h) for h in (Fraction(1, 200), Fraction(9, 1000), Fraction(1, 10000))
            for m, N in ((0, 0), (1, 3), (5, 11), (20, 2))]
    return dict(rows=rows, checks={name: all(r['checks'][name] for r in rows)
                                  for name in rows[0]['checks']})


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != expected_protocol():
        raise ValueError('Frozen scope or parameters changed')
    reviews = dict(all_stage_potential_and_fixed_jet_orders=all_stage_order_review(),
                   actual_q_cutoff_geometry=cutoff_geometry_review(),
                   complete_curl_cutoff_product_rule=curl_review(),
                   complete_cutoff_residual=residual_review(),
                   whole_interval_scale_choice=scale_choice_review(),
                   fixed_stage_flat_residual_comparison=comparison_review(),
                   protected_core_and_local_energy=core_energy_review(),
                   raw_series_and_quantifier_controls=divergent_series_review())
    inherited = inherited_provenance()
    gates = {name: all(reviews[name]['checks'].values()) for name in REVIEW_NAMES}
    gates['all_inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['conditional_local_result_not_promoted'] = not any(FLAGS[k] for k in
        ('actual_global_profile_and_field_numerically_evaluated',
         'actual_infinite_stage_coefficients_or_cutoff_scales_evaluated',
         'inherited_foundation_and_global_auxiliary_geometry_independently_verified',
         'closed_time_auxiliary_atlas_and_one_sided_potential_jets_independently_verified',
         'raw_uncut_infinite_series_convergence_verified',
         'unconditional_complete_residual_construction_verified',
         'admissible_smooth_compact_force_verified', 'whole_space_field_energy_verified',
         'global_solution_exclusion_argument_verified', 'blowup_verified'))
    return encode(dict(status=STATUS, protocol=expected_protocol(), flags=FLAGS,
                       conditional_proof=proof_ledger(), reviews=reviews,
                       inherited_provenance=inherited, source_hashes=source_hashes(),
                       gates=gates, remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS or record.get('conditional_proof') != proof_ledger():
        raise ValueError('Conditional local realization is not an unconditional blow-up construction')
    if record.get('protocol') != expected_protocol() or record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('Frozen scope or remaining obligations changed')
    if set(record.get('gates', {})) != set(GATE_NAMES) or not all(record['gates'].values()):
        raise ValueError('A required cutoff, summation, residual or provenance gate failed')


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
            raise ValueError('Frozen diagnostic differs at '+location)
    elif type(actual) is not type(frozen) or actual != frozen:
        raise ValueError('Frozen exact arithmetic or metadata differs at '+location)


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
    print(str(len(record['gates']))+' gates pass; full-potential realization and local residual flatness remain conditional.')
    print('Endpoint geometry, admissible global forcing, whole-space conclusions and independent foundation review remain open.')


if __name__ == '__main__':
    main()
