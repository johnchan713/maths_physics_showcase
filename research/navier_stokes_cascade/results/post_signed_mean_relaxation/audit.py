#!/usr/bin/env python3
"""Replay post-signed full angular-mean relaxation; wave residual stays open."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from relax import load,fast_review,slow_review,pressure_review,exponent_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='actual-post-signed-fast-slow-angular-mean-relaxation-audited'
FLAGS=dict(current_post_signed_fast_and_slow_mean_inputs_joined=True,
           full_angular_mean_tangential_residual_order_3_minus_6kappa=True,
           recomputed_three_defects_order_3_minus_3kappa=True,
           inherited_nonzero_slow_Haar_mean_and_fast_correlations_retained=True,
           complete_axial_pressure_moment_and_both_base_tails_retained=True,
           actual_flat_velocities_kept_in_every_product=True,
           mean_increment_physical_velocity_pressure_and_partial_curl_potential_bounds_derived=True,
           generated_mean_wave_sources_supported_with_order_2_minus_4kappa=True,
           relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas=True,
           manufactured_exact_checks_not_actual_field_evaluation=True,
           actual_giant_profile_or_first_band_numerically_evaluated=False,
           old_slow_Haar_mean_asserted_zero=False,
           slow_updated_defects_asserted_to_have_order_2alpha=False,
           mean_maps_asserted_to_set_nonlinear_defects_zero=False,
           current_nonzero_angular_residual_cancelled=False,
           every_later_complex_wave_source_and_full_cycle_verified=False,
           unconditional_complete_residual_improvement_verified=False,
           complete_infinite_potential_realization_verified=False,
           blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['audit_current_nonzero_angular_residual_class_support_and_first_sourced_pulse_inverse',
           'verify_later_general_complex_covariance_mean_cycles_without_initial_quadrature_shortcut',
           'join_complete_residual_viscosity_cutoffs_and_all_interactions_at_every_finite_stage',
           'establish_infinite_realization_forcing_and_energy_before_any_blowup_claim']
DEPENDENCIES=('real_signed_covariance_update','complete_covariance_initial_means','compact_mean_corrections',
              'physical_wave_interface','residual_cycle_closure','physical_pulse_transfer',
              'coefficient_support_audit','geometric_pulse_replay')
REVIEWERS=dict(fast_inverse_and_nonzero_old_new_correlations=fast_review,
               slow_recomputed_defects_with_existing_Haar_mean=slow_review,
               all_three_pressures_and_integrated_fluxes=pressure_review,
               all_fixed_order_mean_gain_and_physical_ledger=exponent_review)
EXTRA_GATES={'actual_signed_parent_replayed','frozen_provenance_preserved'}


def protocol():
    return dict(parent_commit='a68e4e5d204e51dcc9ec6d561e10fecaf4985445',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['8.3','8.12-8.17','8.20-8.27'],
                same_Md64_h_N_positive_lambda_width_and_global_geometry=True,
                kappa='1/100000',alpha='2-3*kappa',H='2-4*kappa',nu='3-6*kappa',mu='1-kappa',
                old_Haar_mean_orders=['2 radial','1 swirl','1 axial'],
                final_full_angular_mean_tangential_order='3-6*kappa',final_three_defect_order='3-3*kappa',
                operation_order='current fast inverse; recompute pressure and defects; current slow five-map; recompute again',
                updated_full_six_entry_wave_covariance_held_fixed_during_both_mean_steps=True,
                actual_A0_A1_and_flat_base_residual_retained=True,
                whole_moving_pressure_bump_product_differentiated=True,
                nonzero_old_Haar_mean_not_replaced_by_initial_zero_mean=True,
                partial_radial_axial_mean_potential_not_a_full_potential_realization=True,
                fixture='exact rational one-circle graph and radially compact finite-regularity polynomials',
                fixture_graph_M='3/7',fixture_graph_dr=2,fixture_epsilon='1/200',fixture_radial_interval=[1,2],
                fixture_temporal_source_check_not_actual_global_residual_solve=True,
                actual_two_torus_and_fixed_five_map_estimates_retained_as_analytical_lemmas=True,
                diagnostics_not_continuum_proofs=True)


@lru_cache(None)
def parent_module():
    p=HERE.parent/'real_signed_covariance_update';sys.path.insert(0,str(p))
    try:
        return load('post_signed_parent_real_update',p/'audit.py')
    finally:
        sys.path.pop(0)


@lru_cache(None)
def parent_review():
    parent=parent_module();record=parent.generate();parent.validate(record)
    frozen=json.loads((HERE.parent/'real_signed_covariance_update/evidence.json').read_text())
    parent.compare_records(record,frozen)
    inherited=dict(record['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['real_signed_covariance_update:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return record,inherited


def source_hashes():
    paths=[HERE/n for n in ('README.md','relax.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol(): raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review()
    reviews={n:f() for n,f in REVIEWERS.items()}
    gates={n:bool(r['checks']) and all(r['checks'].values()) for n,r in reviews.items()}
    gates['actual_signed_parent_replayed']=all(parent['gates'].values()) and all(parent['flags'][n] for n in
        ('actual_weighted_fixed_primary_signed_inverse_inputs_joined',
         'first_signed_auxiliary_averaged_tangential_residual_order_3_minus_6kappa',
         'all_six_covariance_entries_and_axial_mean_pressure_recomputed'))
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
                inherited_provenance=dict(inherited),source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):
    return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:
        raise ValueError('Changed construction, existing slow mean or restricted scope')
    if r.get('remaining_obligations')!=REMAINING: raise ValueError('Outstanding obligation changed')
    if set(r.get('reviews',{}))!=set(REVIEWERS) or set(r.get('gates',{}))!=set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(r['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in r['reviews'].values()):
        raise ValueError('Mean relaxation check failed')
    if len(r.get('inherited_provenance',{}))!=625 or not all(r['inherited_provenance'].values()):
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
    print('Angular means only; old nonzero wave residual and later complex cycles/full residual remain open.')


if __name__=='__main__': main()
