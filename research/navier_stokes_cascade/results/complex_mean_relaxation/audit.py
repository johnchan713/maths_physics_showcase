#!/usr/bin/env python3
"""Fresh fast/slow mean maps at the current complex-wave state."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from relaxation import load,fast_review,slow_review,pressure_review,exponent_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-complex-state-full-angular-mean-relaxation-audited'
FLAGS=dict(current_complete_complex_wave_covariance_held_fixed_through_both_mean_maps=True,
    current_fast_then_slow_mean_inputs_joined=True,
    full_angular_mean_tangential_order_2_minus_6kappa_derived=True,
    recomputed_three_defects_order_5_over_2_minus_3kappa_derived=True,
    accumulated_previous_slow_mean_and_old_new_fast_correlations_retained=True,
    postfast_targets_required_by_explicit_stale_order_inequality=True,
    complete_axial_pressure_and_both_background_tails_retained=True,
    actual_flat_velocities_in_every_product_and_Gaussian_tail_exclusion_preserved=True,
    physical_increment_jets_and_added_wave_source_supports_derived=True,
    existing_scoped_replay_adapter_preserved_without_relaxing_any_further_field=True,
    relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas=True,
    manufactured_exact_fixtures_not_actual_global_field_evaluation=True,
    current_Haar_mean_replaced_by_original_initial_value=False,
    old_slow_Haar_mean_asserted_zero=False,
    postslow_defects_asserted_zero_or_order_2alpha=False,
    fast_viscosity_asserted_to_limit_current_full_mean=False,
    old_nonzero_regular_residual_claimed_improved_or_cancelled=False,
    actual_giant_profile_or_two_torus_inverse_numerically_evaluated=False,
    complete_direct_swirl_or_infinite_potential_realization_verified=False,
    all_later_finite_cycles_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['audit_next_current_regular_source_inverse_at_order_3_over_2_minus_6kappa_with_exact_tail_exclusion',
    'join_generated_sources_covariance_and_current_mean_maps_in_a_full_finite_stage_induction',
    'establish_infinite_realization_forcing_and_energy_before_any_blowup_claim']
DEPENDENCIES=('complex_signed_covariance_update','post_signed_mean_relaxation','first_sourced_wave_update',
    'compact_mean_corrections','physical_pulse_transfer','geometric_pulse_replay',
    'residual_cycle_closure','physical_wave_interface')
REVIEWERS=dict(current_fast_inverse_and_accumulated_Haar_state=fast_review,
    current_slow_defects_and_stale_target_obstruction=slow_review,
    complete_fixed_complex_tensor_pressure_and_velocity_moments=pressure_review,
    current_mean_gain_physical_and_uncancelled_wave_ledger=exponent_review)
EXTRA_GATES={'current_generic_signed_parent_replayed','frozen_provenance_preserved'}


def protocol():
    return dict(parent_commit='71965da8e3759b0db0d02b05c298858121c61d10',
        paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_equations=['8.12-8.27','9.3-9.10'],same_Md64_h_N_positive_lambda_width_and_global_geometry=True,
        kappa='1/100000',alpha='3/2-3*kappa',H='3/2-4*kappa',nu='2-6*kappa',mu='1-kappa',
        unchanged_nonzero_order='3/2-6*kappa',final_defect_order='5/2-3*kappa',
        fast_defect_change_order='5/2-5*kappa',stale_order_shortfall='2*kappa',
        old_Haar_mean_orders=['2 radial','1 swirl','1 axial'],
        old_Haar_mean_is_initial_slow_plus_previous_slow_increment=True,
        operation_order='current fast map; recompute pressure and defects; current slow five-map; recompute again',
        complete_current_complex_wave_tensor_fixed=True,
        full_mean_limited_by_old_averaged_covariance_not_fast_viscosity=True,
        all_six_covariances_both_base_tails_actual_A0_A1_and_full_axial_pressure_retained=True,
        designated_Gaussian_flat_tails_never_forward_inverted=True,
        parent_CI_run=38936441478,parent_CI_conclusion='success',
        fixture='fresh exact rational one-circle means after previous slow map; complete post-complex-signed tensor',
        fixture_finite_endpoint_regularity_not_actual_Cinfinity=True,
        physical_sample_h='1/200',enlarged_chart_domain='4Q<qstar',
        partial_radial_axial_curl_potential_not_complete_swirl_potential=True,
        diagnostics_not_continuum_proofs=True)


@lru_cache(None)
def parent_module():
    p=HERE.parent/'complex_signed_covariance_update';sys.path.insert(0,str(p))
    try: return load('complex_mean_current_signed_parent',p/'audit.py')
    finally: sys.path.pop(0)


@lru_cache(None)
def parent_review():
    parent=parent_module();record=parent.generate();parent.validate(record)
    frozen=json.loads((HERE.parent/'complex_signed_covariance_update/evidence.json').read_text())
    parent.validate(frozen);parent.compare_records(record,frozen)
    inherited=dict(record['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['complex_signed_covariance_update:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return record,inherited


def source_hashes():
    paths=[HERE/n for n in ('README.md','relaxation.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol(): raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review();reviews={n:f() for n,f in REVIEWERS.items()}
    gates={n:bool(r['checks']) and all(r['checks'].values()) for n,r in reviews.items()}
    gates['current_generic_signed_parent_replayed']=all(parent['gates'].values()) and all(parent['flags'][n] for n in
        ('current_complex_old_wave_remainder_used_not_initial_real_shortcut',
         'auxiliary_averaged_tangential_order_2_minus_6kappa_derived',
         'full_six_entry_covariance_axial_pressure_and_three_defects_recomputed'))
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'): return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:
        raise ValueError('Changed current mean state, operation order or restricted scope')
    if r.get('remaining_obligations')!=REMAINING: raise ValueError('Outstanding obligation changed')
    if set(r.get('reviews',{}))!=set(REVIEWERS) or set(r.get('gates',{}))!=set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(r['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in r['reviews'].values()):
        raise ValueError('Current complex-state mean check failed')
    if len(r.get('inherited_provenance',{}))!=671 or not all(r['inherited_provenance'].values()):
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
    print('Current mean step only; nonzero regular residual unchanged and later finite cycles remain open.')


if __name__=='__main__': main()
