#!/usr/bin/env python3
"""Join current all-finite foundations to one complete local realization."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from realized_bounds import load,representative_review,comparison_review,scale_review,scope_review
from realized_checks import reviews as exact_reviews

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-source-single-local-realization-relative-to-retained-core-and-finite-analytic-lemmas'
FLAGS=dict(current_all_finite_source_foundation_is_the_seed=True,
    complete_wave_mean_potentials_direct_swirl_and_fixed_gauge_pressure_joined=True,
    one_cutoff_sequence_handles_arbitrary_stage_constants_and_every_fixed_jet=True,
    one_real_smooth_divergence_free_local_field_for_q_positive_supplied_relative_to_retained_lemmas=True,
    complete_residual_flat_in_every_fixed_Cartesian_jet_relative_to_retained_lemmas=True,
    complete_curl_pressure_viscosity_and_all_nonlinear_terms_retained=True,
    protected_inner_field_and_exterior_preserved=True,
    exact_Gaussian_ledger_kept_at_one_fixed_comparison_stage=True,
    historical_frozen_evidence_and_parameters_preserved=True,
    manufactured_exact_finite_checks_are_not_continuum_suprema=True,
    retained_exact_core_and_finite_existence_inverse_gap_lemmas_still_dependencies=True,
    actual_infinite_coefficients_or_cutoff_scales_evaluated=False,
    compact_full_direct_swirl_vector_potential_assumed=False,
    uniform_stage_constants_or_raw_uncut_series_convergence_proved=False,
    arbitrary_j_of_q_substitution_justified=False,
    closed_time_endpoint_atlas_and_one_sided_jets_joined_from_current_foundations=False,
    global_localization_initialization_force_energy_and_comparison_joined=False,
    all_retained_analytic_dependencies_discharged=False,
    unconditional_complete_residual_or_blowup_verified=False,
    formal_certificate=False)
REMAINING=['join_closed_time_auxiliary_atlas_and_complete_one_sided_potential_pressure_jets_to_this_same_current_realization',
    'then_join_localization_initialization_admissible_force_whole_space_energy_and_global_comparison_without_dropping_any_term',
    'discharge_retained_exact_core_and_finite_analytic_existence_inverse_gap_dependencies_before_unconditional_claims']
DEPENDENCIES=('finite_induction_foundation','current_complex_foundation','regular_source_foundation','physical_pulse_foundation',
    'compact_mean_foundation','current_potential_realization','full_potential_realization','background_summation')

def protocol():
    return dict(parent_commit='07b30739332662b1bca55018fa735bb0391ae5d5',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        unchanged_parameters='Md64; unit viscosity; same h,N,positive lambda,width C_setup^-20,roots,geometry,carriers and operation order',
        current_hierarchy='B_j=3/2-2*kappa+j*(1/2-kappa); delta=1/2-kappa; B_(j-1)=delta*j+rho; rho=1-kappa',
        representatives='Delta u_j=curl(A_j)+b_j*e_theta; partial_theta b_j=0; full p_j=p^[j]-p^[j-1] in fixed compact gauge',
        physical_losses='ell_m=2*A+2*m; L_m=ell_(m+1); H_m=L_m+1; K_m=2*A+3/2+2*m',
        diagonal='one a_j>=2*a_(j-1); finite lists through m=j; whole interval 0<q<=1/a_j; g_j=h*delta*j',
        tail='2^(-J)*q^(h*delta*(J+1)/2-L_m) for q<=1/(2*a_J), fixed J>=max(1,m)',
        comparison='fix m,N then J>=m+2 with h*delta*(J+1)/2>=N+1+L_(m+2)+H_(m+2) and h*B_J-K_m>=N+1',
        full_residual='e_t-Delta e+grad pi+u_J.grad e+e.grad u_J+e.grad e',
        pressure='full telescoping differences after all four operations, both background Pi1/Pi2, axial memory; gauge fixed before cutoff',
        domain='same 0<q<qstar<1, t<1; no one-sided endpoint or global force conclusion',
        scope='single local realization and complete all-fixed-jet flat residual relative to retained exact-core and finite analytic lemmas')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'finite_induction_foundation';sys.path.insert(0,str(p))
    try:return load('realization_foundation_finite_parent',p/'audit.py')
    finally:sys.path.pop(0)

@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'finite_induction_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,d in f['source_hashes'].items():history['finite_induction_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    return r,history

def sources():
    paths=[HERE/n for n in ('README.md','realized_bounds.py','realized_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'finite_induction_foundation/finite_bounds.py',HERE.parent/'background_summation/summation.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():return deepcopy(dict(current_complete_representatives=representative_review(),
    all_interval_diagonal=scale_review(),fixed_stage_complete_residual_comparison=comparison_review(),
    local_field_support_pressure_and_quantifier_scope=scope_review(),**exact_reviews()))

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews();gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(current_finite_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Current local realization scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Endpoint/global/analytic dependencies erased')
    if set(r.get('gates',{}))!=set(reviews())|{'current_finite_parent_reproduced','historical_bytes_preserved'} or not all(r['gates'].values()):raise ValueError('Realization gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=961 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Current local realization relative to retained lemmas; endpoint, admissible global force and blow-up remain open.')

if __name__=='__main__':main()
