#!/usr/bin/env python3
"""Replay the actual-primary INITIAL mean update; later cycles stay open."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from means import load, covariance_review, haar_update_review, exponent_review, pressure_tail_review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'actual-primary-complete-covariance-initial-mean-gain-audited'
FLAGS = dict(actual_primary_full_six_entry_covariance_initial_mean_inputs_joined=True,
             initial_fast_update_has_exact_zero_Haar_velocity=True,
             initial_recomputed_three_compatibility_defects_in_S2=True,
             initial_angular_mean_annular_residual_in_M2_minus_3kappa_modulo_flat=True,
             complete_axial_velocity_and_moving_pressure_moment_retained=True,
             both_background_pressure_tails_retained_via_same_base_subtraction=True,
             relies_on_replayed_background_geometry_pulse_and_compact_mean_lemmas=True,
             manufactured_rational_checks_not_actual_profile_evaluation=True,
             actual_giant_profile_or_first_band_numerically_evaluated=False,
             subsequent_signed_covariance_waves_and_full_mean_hierarchy_verified=False,
             all_nonlinear_generated_source_class_estimates_verified=False,
             unconditional_complete_residual_improvement_verified=False,
             blowup_verified=False, formal_proof_assistant_certificate=False)
REMAINING = ['realize_signed_mean_covariance_targets_and_bound_complete_updated_curl_covariance',
             'verify_every_nonlinear_generated_source_class_envelope_and_endpoint_extension',
             'close_every_finite_mean_wave_cycle_with_actual_axial_pressure_and_all_interactions',
             'prove_full_residual_viscosity_cutoff_and_physical_scale_bounds_before_any_blowup_claim']
DEPENDENCIES = ('physical_pulse_transfer', 'compact_mean_corrections', 'coefficient_support_audit',
                'residual_cycle_closure', 'physical_wave_interface', 'weighted_background_edges',
                'background_jet_foundation', 'geometric_pulse_replay')
REVIEWERS = dict(exact_complete_curl_and_pressure_moments=covariance_review,
                 exact_Haar_and_recomputed_defect_changes=haar_update_review,
                 all_fixed_order_initial_exponent_ledger=exponent_review,
                 complete_background_pressure_tail_controls=pressure_tail_review)
EXTRA_GATES = {'actual_primary_parent_replayed', 'retained_compact_mean_maps_replayed',
               'frozen_provenance_preserved'}


def protocol():
    return dict(parent_commit='0cb55230b7b80332266cab01daed83d3e336c79e',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['8.3', '8.12-8.17', '8.20-8.27'],
                same_Md64_h_N_positive_lambda_width_and_global_geometry=True,
                initial_wave='fixed actual complete real-primary curl; no signed covariance update yet',
                order='pressure; fast means; recompute pressure and defects; slow five-map; recompute again',
                kappa='1/100000', raw_defect_order='1', final_three_defect_order='2',
                final_tangential_mean_residual_order='2-3*kappa',
                retained_flat_base_residual_separate_from_weighted_shell_class=True,
                actual_flat_velocity_kept_in_all_products=True,
                first_two_background_pressure_tails='retained; cancel only in SAME full base subtraction',
                auxiliary_average_not_angular_average=True,
                every_fixed_coefficient_derivative_order=True,
                initial_zero_Haar_cancellation_not_a_later_cycle_estimate=True,
                fixture_arithmetic='Fraction rational polynomial and degree-two Haar moments',
                fixture_radial_support=[1, 2], fixture_epsilon='1/200', fixture_angular_carrier=3,
                finite_fixture_endpoint_regularity_not_actual_Cinfinity_support=True,
                diagnostics_are_not_continuum_proofs=True)


@lru_cache(None)
def module(folder):
    path = HERE.parent/folder
    sys.path.insert(0, str(path))
    try:
        return load('initial_mean_dependency_'+folder, path/'audit.py')
    finally:
        sys.path.pop(0)


@lru_cache(None)
def parent_review():
    parent = module('physical_pulse_transfer')
    record = parent.generate(); parent.validate(record)
    frozen = json.loads((HERE.parent/'physical_pulse_transfer/evidence.json').read_text())
    parent.compare_records(record, frozen)
    inherited = dict(record['inherited_provenance'])
    for path, digest in frozen['source_hashes'].items():
        inherited['physical_pulse_transfer:'+path] = hashlib.sha256((PROJECT/path).read_bytes()).hexdigest() == digest
    return record, inherited


@lru_cache(None)
def compact_review():
    compact = module('compact_mean_corrections')
    record = compact.generate(); compact.validate(record)
    frozen = json.loads((HERE.parent/'compact_mean_corrections/evidence.json').read_text())
    compact.compare_records(record, frozen)
    return dict(status=record['status'], gates=record['gates'],
                inherited_conditional_scope_not_silently_changed=True)


def source_hashes():
    paths = [HERE/n for n in ('README.md', 'means.py', 'audit.py', 'test_audit.py',
                              'protocol.json', 'requirements.txt', '.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != protocol():
        raise ValueError('Frozen protocol changed')
    parent, inherited = parent_review()
    compact = compact_review()
    reviews = {name: reviewer() for name, reviewer in REVIEWERS.items()}
    gates = {name: bool(r['checks']) and all(r['checks'].values()) for name, r in reviews.items()}
    gates['actual_primary_parent_replayed'] = all(parent['gates'].values()) and parent['flags'][
        'actual_background_inputs_to_primary_pulse_theorem_discharged']
    gates['retained_compact_mean_maps_replayed'] = all(compact['gates'].values())
    gates['frozen_provenance_preserved'] = bool(inherited) and all(inherited.values())
    return dict(status=STATUS, protocol=protocol(), flags=dict(FLAGS), reviews=reviews,
                compact_map_replay=compact, gates=gates, inherited_provenance=dict(inherited),
                source_hashes=source_hashes(), remaining_obligations=REMAINING[:])


def compare_records(actual, frozen, path='root'):
    return module('physical_pulse_transfer').compare_records(actual, frozen, path)


def validate(record):
    if record.get('status') != STATUS or record.get('protocol') != protocol() or record.get('flags') != FLAGS:
        raise ValueError('Changed construction, initial-only scope or analytical dependencies')
    if record.get('remaining_obligations') != REMAINING:
        raise ValueError('Outstanding obligation changed')
    if set(record.get('reviews', {})) != set(REVIEWERS) or set(record.get('gates', {})) != set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(record['gates'].values()) or not all(r.get('checks') and all(r['checks'].values()) for r in record['reviews'].values()):
        raise ValueError('Initial mean check failed')
    if len(record.get('inherited_provenance', {})) != 595 or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited provenance changed')
    if record.get('source_hashes') != source_hashes():
        raise ValueError('Current sources changed')
    compare_records(parent_review()[1], record['inherited_provenance'], 'inherited_provenance')
    compare_records(compact_review(), record.get('compact_map_replay'), 'compact_map_replay')
    for name, reviewer in REVIEWERS.items():
        compare_records(reviewer(), record['reviews'][name], name)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify-record', type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit('Output exists; choose a fresh path')
    record = generate(); validate(record)
    if args.verify_record:
        frozen = json.loads(args.verify_record.read_text()); validate(frozen); compare_records(record, frozen)
    args.output.write_text(json.dumps(record, sort_keys=True, indent=2)+'\n')
    print(STATUS)
    print(f"{len(record['gates'])} gates pass; {len(record['inherited_provenance'])} inherited hashes agree.")
    print('Initial means only; signed covariance waves, later source closure, full residual and blow-up remain open.')


if __name__ == '__main__':
    main()
