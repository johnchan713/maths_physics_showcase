#!/usr/bin/env python3
"""Current-seeded conditional complete finite-stage induction milestone."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from induction import load,seed_review,half_line_review,physical_quantifier_review,mechanism_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='conditional-current-state-complete-finite-stage-induction'
FLAGS=dict(current_audited_complex_mean_state_supplies_actual_seed=True,
    conditional_complete_four_operation_finite_stage_induction_derived=True,
    positive_cycle_order_gain_one_half_minus_3kappa_on_whole_half_line=True,
    exact_gain_limit_from_full_linear_and_old_complex_cross_derived=True,
    every_generated_regular_source_support_Gaussian_envelope_and_endpoint_hypothesis_preserved=True,
    actual_flat_velocities_and_exact_excluded_Gaussian_tail_ledger_retained=True,
    all_six_covariances_full_axial_pressure_both_base_tails_and_current_defects_retained=True,
    reality_divergence_and_two_Haar_velocity_moments_preserved=True,
    accumulated_old_slow_Haar_mean_order_one_preserved=True,
    same_fixed_inverse_domain_and_stage_independent_physical_jet_loss_derived_conditionally=True,
    for_each_fixed_physical_power_and_jet_order_a_finite_state_supplied_conditionally=True,
    finite_stage_induction_stream_reaches_conditional_stopping_milestone=True,
    relies_on_replayed_leading_background_geometry_and_inverse_analytical_premises=True,
    inherited_manufactured_controls_not_continuum_proofs=True,
    independent_foundation_review_complete=False,
    actual_giant_fields_or_thousand_stage_sequence_numerically_evaluated=False,
    initial_real_quadrature_shortcut_used_for_generic_update=False,
    designated_Gaussian_tails_forward_inverted_or_erased=False,
    uniform_stage_constants_or_fixed_q_norm_contraction_verified=False,
    arbitrary_j_of_q_substitution_justified=False,
    single_infinite_realized_field_with_flat_complete_residual_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    admissible_global_smooth_force_and_complete_energy_verified=False,
    blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['join_complete_potential_representatives_and_whole_interval_diagonal_realization_to_current_finite_states',
    'verify_one_sided_endpoint_localization_admissible_smooth_force_energy_and_surviving_lower_bound',
    'independently_review_inherited_leading_background_and_global_auxiliary_geometry_before_unconditional_claims']
DEPENDENCIES=('complex_mean_relaxation','complex_signed_covariance_update','first_sourced_wave_update',
    'post_signed_mean_relaxation','actual_band_pulse_bounds','compact_mean_corrections',
    'residual_cycle_closure','physical_pulse_transfer')


def protocol():
    return dict(parent_commit='e3b5dd1fd65003b75bc33181bad0711b40a1c9b5',
        paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_equations=['6.19-6.32','7.39-7.42','8.12-8.27','9.1-9.10','Propositions 9.3-9.6'],
        same_Md64_h_N_positive_lambda_width_and_global_geometry=True,kappa='1/100000',
        initial_B='3/2-6*kappa',gain='1/2-3*kappa',rho='1-3*kappa',mu='1-kappa',
        state_mean_and_defect_order='B+1/2',postslow_defect_order='B+3/2',
        cycle_order='regular sourced wave; current generic signed covariance; current fast mean; current slow five-map',
        pressure_and_defects_recomputed_after_every_velocity_change=True,
        designated_Gaussian_tails_never_forward_inverted=True,
        exact_flat_velocities_kept_in_all_actual_products=True,
        proof_domain='all B>=3/2-6*kappa by affine slopes and boundary margins, not a finite stage scan',
        physical_loss='K_m=2*A+3/2+2*m',stage_order='B_j=3/2-6*kappa+j*(1/2-3*kappa)',
        harmonic_bound='4*2^j',common_inverse_domain='same qstar; enlarged charts 4Q<qstar',
        constants_may_grow_arbitrarily_in_finite_stage=True,
        parent_CI_run=38944616973,parent_CI_conclusion='success',
        physical_sample_h='1/200',
        quantifier='for each fixed target and finite jet order exists a finite state; no single-field conclusion',
        diagnostics_not_continuum_proofs=True)


@lru_cache(None)
def parent_module():
    p=HERE.parent/'complex_mean_relaxation';sys.path.insert(0,str(p))
    try: return load('finite_induction_current_mean_parent',p/'audit.py')
    finally: sys.path.pop(0)


@lru_cache(None)
def parent_review():
    parent=parent_module();record=parent.generate();parent.validate(record)
    frozen=json.loads((HERE.parent/'complex_mean_relaxation/evidence.json').read_text())
    parent.validate(frozen);parent.compare_records(record,frozen)
    inherited=dict(record['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['complex_mean_relaxation:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return record,inherited


REVIEWERS=dict(current_actual_seed=lambda:seed_review(parent_review()[0]),
    complete_four_step_half_line_induction=half_line_review,
    stage_independent_physical_loss_and_quantifier_boundary=physical_quantifier_review,
    replayed_complete_interaction_and_missing_hypothesis_controls=mechanism_review)
EXTRA_GATES={'current_complex_mean_parent_replayed','frozen_provenance_preserved'}


def source_hashes():
    paths=[HERE/n for n in ('README.md','induction.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol(): raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review();reviews={n:f() for n,f in REVIEWERS.items()}
    gates={n:bool(r['checks']) and all(r['checks'].values()) for n,r in reviews.items()}
    gates['current_complex_mean_parent_replayed']=all(parent['gates'].values()) and all(parent['flags'][n] for n in
        ('current_fast_then_slow_mean_inputs_joined','full_angular_mean_tangential_order_2_minus_6kappa_derived',
         'recomputed_three_defects_order_5_over_2_minus_3kappa_derived',
         'actual_flat_velocities_in_every_product_and_Gaussian_tail_exclusion_preserved'))
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'): return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:
        raise ValueError('Changed conditional induction, seed, tail ledger or quantifier scope')
    if r.get('remaining_obligations')!=REMAINING: raise ValueError('Outstanding obligation changed')
    if set(r.get('reviews',{}))!=set(REVIEWERS) or set(r.get('gates',{}))!=set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(r['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in r['reviews'].values()):
        raise ValueError('Conditional finite induction check failed')
    if len(r.get('inherited_provenance',{}))!=686 or not all(r['inherited_provenance'].values()):
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
    print('Conditional finite-stage milestone only; single-field realization, force and blow-up remain open.')


if __name__=='__main__': main()
