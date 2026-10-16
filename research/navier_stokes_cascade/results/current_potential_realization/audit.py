#!/usr/bin/env python3
"""Join the current finite-state hierarchy to one conditional local realization."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from realization import load,order_review,comparison_review,scale_review,identity_review,scope_review,old_module

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='conditional-current-state-complete-potential-diagonal-realization'
FLAGS=dict(current_audited_finite_state_hierarchy_is_the_seed=True,
    complete_wave_mean_curl_direct_swirl_and_pressure_representatives_joined=True,
    one_scale_sequence_handles_arbitrary_finite_stage_constants_and_all_fixed_jets=True,
    conditional_single_real_smooth_divergence_free_local_field_for_q_positive_supplied=True,
    conditional_complete_residual_flatness_in_every_fixed_Cartesian_jet_supplied=True,
    exact_protected_inner_field_and_exterior_preserved_conditionally=True,
    designated_Gaussian_tails_never_inverted_or_independently_summed=True,
    historical_frozen_evidence_preserved_and_replayed=True,
    inherited_manufactured_diagnostics_are_not_continuum_proofs=True,
    relies_on_inherited_leading_background_geometry_and_inverse_analytical_premises=True,
    actual_infinite_sequence_coefficients_or_scales_evaluated=False,
    direct_swirl_has_a_proved_compact_full_vector_potential=False,
    uniform_stage_constants_raw_uncut_series_or_fixed_q_contraction_verified=False,
    independent_foundation_review_complete=False,
    closed_time_endpoint_atlas_and_one_sided_jets_independently_verified=False,
    physical_localization_initialization_admissible_global_force_and_energy_verified=False,
    unconditional_complete_residual_construction_verified=False,
    blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['audit_closed_time_auxiliary_atlas_and_complete_one_sided_potential_pressure_jets_away_from_singular_corner',
    'join_complete_physical_localization_and_initialization_to_current_realization_with_every_curl_pressure_viscous_and_nonlinear_term',
    'verify_admissible_smooth_compact_force_whole_space_energy_and_global_solution_exclusion',
    'independently_review_inherited_leading_background_inverse_and_global_geometry_before_unconditional_claims']
DEPENDENCIES=('finite_stage_induction','full_potential_realization','complex_mean_relaxation',
    'first_sourced_wave_update','physical_pulse_transfer','compact_mean_corrections',
    'background_summation','actual_band_pulse_bounds')


def protocol():
    return dict(parent_commit='d065d09fa2f16710ea590454ebe5389e5edfb909',
        paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_equations=['9.15-9.21','Lemma 9.8','Proposition 9.9'],
        actual_finite_hierarchy='B_j=3/2-6*kappa+j*(1/2-3*kappa)',
        increment_gain='g_j=h*(1/2-3*kappa)*j',tail_gain='g_j/2',
        representative='Delta u_j=curl(A_j)+b_j*e_theta; partial_theta b_j=0',
        complete_pressure='p_j=p^[j]-p^[j-1] after all four operations',
        uniform_physical_losses='ell_m=2*A+2*m; L_m=ell_(m+1); H_m=L_m+1; K_m=2*A+3/2+2*m',
        scale_requirements='finite list of physical jets through m=j; whole interval 0<q<=1/a_j',
        fixed_comparison='J>=m+2; h*delta*(J+1)/2>=N+1+L_(m+2)+H_(m+2); h*B_J-K_m>=N+1',
        domain='same inherited 0<q<qstar<1 and t<1; no closed-time endpoint assertion',
        same_h_N_lambda_width_labels_and_fixed_primary_roots=True,
        stage_constants_may_grow_arbitrarily_fast=True,
        flat_tail_ledger_compared_at_one_fixed_stage_never_forward_inverted=True,
        manufactured_diagnostics_not_actual_global_field_evaluation=True,
        claim='conditional local realization and all-fixed-jet complete residual flatness only')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'finite_stage_induction';sys.path.insert(0,str(p))
    try: return load('current_realization_finite_induction_parent',p/'audit.py')
    finally: sys.path.pop(0)


@lru_cache(None)
def parent_review():
    p=parent_module();r=p.generate();p.validate(r)
    frozen=json.loads((HERE.parent/'finite_stage_induction/evidence.json').read_text())
    p.validate(frozen);p.compare_records(r,frozen)
    inherited=dict(r['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['finite_stage_induction:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return r,inherited


@lru_cache(None)
def historical_review():
    p=old_module();r=p.generate();p.validate(r)
    frozen=json.loads((HERE.parent/'full_potential_realization/evidence.json').read_text())
    p.validate(frozen);p.compare_records(r,frozen)
    return dict(old_conditional_hierarchy_not_substituted_for_current=True,
        checks=dict(historical_frozen_record_reproduces=all(r['gates'].values()),
            historical_sources_still_match=r['source_hashes']==frozen['source_hashes']))


REVIEWERS=dict(current_all_stage_representatives_and_physical_losses=order_review,
    whole_interval_diagonal_with_current_gain=scale_review,
    fixed_finite_stage_complete_residual_comparison=comparison_review,
    complete_implicit_q_curl_and_residual_controls=identity_review,
    exact_tail_pressure_and_quantifier_scope=scope_review,
    historical_realization_replayed_not_assumed_current=historical_review)
EXTRA_GATES={'current_actual_finite_hierarchy_replayed','frozen_provenance_preserved'}


def source_hashes():
    paths=[HERE/n for n in ('README.md','realization.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol(): raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review();reviews={k:f() for k,f in REVIEWERS.items()}
    gates={k:bool(r['checks']) and all(r['checks'].values()) for k,r in reviews.items()}
    gates['current_actual_finite_hierarchy_replayed']=all(parent['gates'].values()) and all(parent['flags'][k] for k in
        ('conditional_complete_four_operation_finite_stage_induction_derived',
         'same_fixed_inverse_domain_and_stage_independent_physical_jet_loss_derived_conditionally',
         'actual_flat_velocities_and_exact_excluded_Gaussian_tail_ledger_retained'))
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
        inherited_provenance=inherited,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'): return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:
        raise ValueError('Conditional local realization, seed, representatives or scope changed')
    if r.get('remaining_obligations')!=REMAINING: raise ValueError('Outstanding obligation changed')
    if set(r.get('reviews',{}))!=set(REVIEWERS) or set(r.get('gates',{}))!=set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(r['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in r['reviews'].values()):
        raise ValueError('Conditional realization check failed')
    if len(r.get('inherited_provenance',{}))!=701 or not all(r['inherited_provenance'].values()):
        raise ValueError('Inherited provenance changed')
    if r.get('source_hashes')!=source_hashes(): raise ValueError('Current sources changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')
    for k,f in REVIEWERS.items(): compare_records(f(),r['reviews'][k],'root.reviews.'+k)


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
    print('Conditional local single-field milestone only; endpoint, global force and blow-up remain open.')


if __name__=='__main__': main()
