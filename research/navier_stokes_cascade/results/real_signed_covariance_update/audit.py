#!/usr/bin/env python3
"""Replay first actual phase-locked signed covariance update; not a full cycle."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from signed import load, covariance_review, inverse_review, pressure_review, exponent_review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'actual-first-phase-locked-real-signed-covariance-update-audited'
FLAGS = dict(actual_weighted_fixed_primary_signed_inverse_inputs_joined=True,
             first_phase_locked_real_signed_full_curl_covariance_error_order_3_minus_6kappa=True,
             first_signed_auxiliary_averaged_tangential_residual_order_3_minus_6kappa=True,
             all_six_covariance_entries_and_axial_mean_pressure_recomputed=True,
             this_substep_generated_wave_source_class_and_endpoint_inputs_derived=True,
             signed_velocity_pressure_potential_fixed_physical_jets_transferred=True,
             both_background_pressure_tails_and_actual_flat_velocities_retained=True,
             relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas=True,
             exact_manufactured_diagnostics_not_actual_profile_evaluation=True,
             actual_giant_profile_or_first_band_numerically_evaluated=False,
             sharper_covariance_formula_applies_to_arbitrary_later_complex_amplitudes=False,
             signed_linear_inverse_is_an_exact_nonlinear_covariance_solve=False,
             compatibility_defects_asserted_to_remain_in_S2_after_wave_change=False,
             complete_unaveraged_mean_residual_already_improved=False,
             every_later_nonlinear_source_and_complete_cycle_verified=False,
             unconditional_complete_residual_improvement_verified=False,
             blowup_verified=False, formal_proof_assistant_certificate=False)
REMAINING = ['apply_fast_then_slow_mean_maps_to_current_recomputed_residual_and_defects',
             'audit_later_sourced_wave_and_general_complex_signed_updates_without_initial_quadrature_shortcut',
             'verify_every_finite_cycle_source_support_pressure_viscosity_and_cutoff_bound',
             'establish_full_residual_infinite_realization_forcing_and_energy_before_any_blowup_claim']
DEPENDENCIES = ('complete_covariance_initial_means','physical_pulse_transfer','actual_band_pulse_bounds',
                'physical_wave_interface','compact_mean_corrections','residual_cycle_closure',
                'geometric_pulse_replay','weighted_background_edges')
REVIEWERS = dict(exact_phase_locked_full_curl_covariance=covariance_review,
                 signed_inverse_edge_and_partition_controls=inverse_review,
                 recomputed_axial_pressure_and_fluxes=pressure_review,
                 initial_covariance_source_and_physical_orders=exponent_review)
EXTRA_GATES = {'actual_initial_means_parent_replayed','frozen_provenance_preserved'}


def protocol():
    return dict(parent_commit='e06495cda01638fe17e20f21f5032c606c806099',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['7.31-7.42','8.3','8.12-8.18'],
                same_Md64_h_N_positive_lambda_width_and_global_geometry=True,
                kappa='1/100000', initial_target_order='2-3*kappa',
                signed_velocity_order='3/2-3*kappa', first_real_pair_error_order='3-6*kappa',
                full_tangential_mean_conservative_order='2-4*kappa', recomputed_defect_order='2-3*kappa',
                inverse_denominators='fixed positive primary amplitudes; never updated roots',
                signed_increment_square_and_all_curl_commutators_retained=True,
                real_means_phase_locked_homogeneous_coefficients_not_merely_real_physical_velocity=True,
                no_new_sourced_wave_has_preceded_this_signed_update=True,
                same_phase_in_each_old_new_label=True,
                Haar_average_removes_only_the_subsequent_radial_graph_operator=True,
                fixed_mean_velocity_during_this_wave_substep=True,
                pressure_and_all_three_defects_recomputed_after_wave_change=True,
                actual_flat_velocities_kept_in_all_products=True,
                background_flat_error_separate_from_shell_weight=True,
                physical_loss_space='1/2+h/2',physical_loss_time='1+3h/2',
                all_fixed_jets_on_every_accepted_later_scale=True,
                fixture='exact rational complete-curl polynomials; finite endpoint regularity',
                fixture_radial_interval=[1,2],fixture_frequencies=[3,5],fixture_epsilon='1/200',
                numerical_or_rational_diagnostics_not_continuum_proofs=True)


@lru_cache(None)
def parent_module():
    path = HERE.parent/'complete_covariance_initial_means'
    sys.path.insert(0,str(path))
    try:
        return load('real_signed_parent_initial_means',path/'audit.py')
    finally:
        sys.path.pop(0)


@lru_cache(None)
def parent_review():
    parent = parent_module(); record = parent.generate(); parent.validate(record)
    frozen = json.loads((HERE.parent/'complete_covariance_initial_means/evidence.json').read_text())
    parent.compare_records(record,frozen)
    inherited = dict(record['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['complete_covariance_initial_means:'+path] = hashlib.sha256((PROJECT/path).read_bytes()).hexdigest() == digest
    return record,inherited


def source_hashes():
    paths = [HERE/n for n in ('README.md','signed.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != protocol():
        raise ValueError('Frozen protocol changed')
    parent,inherited = parent_review()
    reviews = {name:reviewer() for name,reviewer in REVIEWERS.items()}
    gates = {name:bool(r['checks']) and all(r['checks'].values()) for name,r in reviews.items()}
    gates['actual_initial_means_parent_replayed'] = all(parent['gates'].values()) and all(parent['flags'][name] for name in
        ('actual_primary_full_six_entry_covariance_initial_mean_inputs_joined',
         'initial_recomputed_three_compatibility_defects_in_S2',
         'initial_angular_mean_annular_residual_in_M2_minus_3kappa_modulo_flat'))
    gates['frozen_provenance_preserved'] = bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
                inherited_provenance=dict(inherited),source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(actual,frozen,path='root'):
    return parent_module().compare_records(actual,frozen,path)


def validate(record):
    if record.get('status') != STATUS or record.get('protocol') != protocol() or record.get('flags') != FLAGS:
        raise ValueError('Changed phase locking, construction or restricted scope')
    if record.get('remaining_obligations') != REMAINING:
        raise ValueError('Outstanding obligation changed')
    if set(record.get('reviews',{})) != set(REVIEWERS) or set(record.get('gates',{})) != set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(record['gates'].values()) or not all(r.get('checks') and all(r['checks'].values()) for r in record['reviews'].values()):
        raise ValueError('Signed covariance check failed')
    if len(record.get('inherited_provenance',{})) != 610 or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited provenance changed')
    if record.get('source_hashes') != source_hashes():
        raise ValueError('Current sources changed')
    compare_records(parent_review()[1],record['inherited_provenance'],'inherited_provenance')
    for name,reviewer in REVIEWERS.items():
        compare_records(reviewer(),record['reviews'][name],name)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit('Output exists; choose a fresh path')
    record = generate(); validate(record)
    if args.verify_record:
        frozen = json.loads(args.verify_record.read_text()); validate(frozen); compare_records(record,frozen)
    args.output.write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
    print(STATUS)
    print(f"{len(record['gates'])} gates pass; {len(record['inherited_provenance'])} inherited hashes agree.")
    print('First phase-locked update only; current means need correction and later complex sources/full residual remain open.')


if __name__ == '__main__':
    main()
