#!/usr/bin/env python3
"""Current complex-wave signed covariance substep and scoped replay repair."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from replay import parent_module
from covariance import tensor_review,pressure_review,exponent_review,replay_review,signed

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-complex-wave-signed-covariance-substep-audited'
FLAGS=dict(current_complex_old_wave_remainder_used_not_initial_real_shortcut=True,
    actual_fixed_primary_signed_inverse_inputs_and_edge_weights_joined=True,
    generic_old_error_cross_new_curl_cross_and_signed_square_retained=True,
    auxiliary_averaged_tangential_order_2_minus_6kappa_derived=True,
    full_six_entry_covariance_axial_pressure_and_three_defects_recomputed=True,
    nonzero_regular_residual_preserved_at_3_over_2_minus_6kappa=True,
    physical_jets_and_generated_source_supports_of_this_substep_joined=True,
    exact_Gaussian_tails_and_actual_flat_velocities_retained=True,
    observed_parent_CI_replay_failure_corrected_by_explicit_three_field_adapter=True,
    all_historical_source_and_evidence_bytes_preserved=True,
    relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas=True,
    manufactured_diagnostics_not_actual_giant_band_evaluation=True,
    parent_CI_claimed_successful=False,
    initial_real_quadrature_shortcut_used_for_generic_update=False,
    all_covariance_entries_claimed_to_gain_to_2_minus_6kappa=False,
    full_means_and_defects_claimed_relaxed_in_this_substep=False,
    old_nonzero_regular_residual_claimed_further_improved=False,
    Gaussian_flat_tails_forward_inverted_or_erased=False,
    actual_giant_profile_or_first_band_numerically_evaluated=False,
    every_later_complex_wave_source_and_full_cycle_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    complete_infinite_potential_realization_verified=False,
    blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['relax_current_full_means_with_fast_then_slow_maps_recomputing_pressure_and_defects_each_time',
    'verify_all_later_sources_supports_and_complete_finite_cycles_with_exact_tail_ledger',
    'establish_infinite_realization_forcing_and_energy_before_any_blowup_claim']
DEPENDENCIES=('first_sourced_wave_update','post_signed_mean_relaxation','real_signed_covariance_update',
    'physical_pulse_transfer','actual_band_pulse_bounds','residual_cycle_closure',
    'compact_mean_corrections','physical_wave_interface')
REVIEWERS=dict(generic_complete_tensor_and_complex_old_error=tensor_review,
    recomputed_full_axial_pressure_and_defects=pressure_review,
    current_fixed_signed_inverse=signed.inverse_review,
    current_generic_orders_and_physical_sources=exponent_review,
    parent_CI_roundoff_compatibility_not_scope_relaxation=replay_review)
EXTRA_GATES={'actual_first_sourced_parent_replayed','frozen_provenance_preserved'}


def protocol():
    return dict(parent_commit='5e7f1ee0103be6f9ba9ad6b1f9e4fec64efffb31',
        paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_equations=['7.39-7.42','8.12-8.27','9.3-9.6'],
        same_Md64_h_N_positive_lambda_width_and_global_geometry=True,
        kappa='1/100000',B='1-3*kappa',alpha='3/2-3*kappa',rho='1-3*kappa',
        nu='2-6*kappa',unchanged_nonzero_order='3/2-6*kappa',
        full_mean_order='3/2-4*kappa',recomputed_defect_order='3/2-3*kappa',
        operation_order='current Haar target; original fixed signed inverse; full curl; recompute six covariances and pressure',
        old_error_cross_and_signed_square_both_limit_generic_gain=True,
        designated_Gaussian_flat_tails_never_forward_inverted=True,
        whole_moving_pressure_bump_product_differentiated=True,
        all_six_covariances_both_base_tails_and_actual_A0_A1_retained=True,
        old_slow_Haar_mean_not_replaced_by_initial_zero=True,
        parent_CI_run=38929673498,parent_CI_conclusion='failure',
        compatibility_adapter_fields=['maximum_relative_commutator_error','maximum_curl_product_error','maximum_divergence_error'],
        compatibility_adapter_namespace='full_viscous_potential_pressure_cutoff',
        compatibility_adapter_absolute_ceiling='2e-12',original_acceptance_gates_unchanged=True,
        fixture='manufactured phase-rotated exact rational complete curls and nonconstant signed potential',
        physical_sample_h='1/200',enlarged_chart_domain='4Q<qstar',
        diagnostics_not_continuum_proofs=True)


@lru_cache(None)
def parent_review():
    parent=parent_module();record=parent.generate();parent.validate(record)
    frozen=json.loads((HERE.parent/'first_sourced_wave_update/evidence.json').read_text())
    parent.validate(frozen);parent.compare_records(record,frozen)
    inherited=dict(record['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['first_sourced_wave_update:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return record,inherited


def source_hashes():
    paths=[HERE/n for n in ('README.md','covariance.py','replay.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol(): raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review();reviews={n:f() for n,f in REVIEWERS.items()}
    gates={n:bool(r['checks']) and all(r['checks'].values()) for n,r in reviews.items()}
    gates['actual_first_sourced_parent_replayed']=all(parent['gates'].values()) and all(parent['flags'][n] for n in
        ('current_regular_source_class_support_and_endpoint_hypotheses_joined',
         'all_six_complex_covariances_and_full_axial_pressure_recomputed',
         'mean_residual_and_three_defects_regenerated_not_preserved'))
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'): return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:
        raise ValueError('Changed current complex construction, replay repair or restricted scope')
    if r.get('remaining_obligations')!=REMAINING: raise ValueError('Outstanding obligation changed')
    if set(r.get('reviews',{}))!=set(REVIEWERS) or set(r.get('gates',{}))!=set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(r['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in r['reviews'].values()):
        raise ValueError('Generic signed covariance check failed')
    if len(r.get('inherited_provenance',{}))!=655 or not all(r['inherited_provenance'].values()):
        raise ValueError('Inherited provenance changed')
    if r.get('source_hashes')!=source_hashes(): raise ValueError('Current sources changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')
    for n,f in REVIEWERS.items(): compare_records(f(),r['reviews'][n],n)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args(argv)
    if args.output.exists(): raise SystemExit('Output exists; choose a fresh path')
    r=generate();validate(r)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text());validate(frozen);compare_records(r,frozen)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(STATUS)
    print(f"{len(r['gates'])} gates pass; {len(r['inherited_provenance'])} inherited hashes agree.")
    print('One generic signed substep; full means, later cycles and complete residual remain open.')


if __name__=='__main__': main()
