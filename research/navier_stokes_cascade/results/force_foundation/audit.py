#!/usr/bin/env python3
"""Join the same current field to compact smooth forcing and preterminal energy."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from localized_bounds import load,domain_review,terminal_review,extension_review,energy_review,scope_review
from localized_checks import reviews as exact_reviews

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-source-compact-smooth-force-and-preterminal-energy-relative-to-retained-core-and-finite-analytic-lemmas'
FLAGS=dict(same_current_single_realization_is_localized=True,
    enlarged_one_quarter_to_four_chart_margin_retained=True,
    axisymmetric_potential_cutoff_and_direct_swirl_preserve_divergence=True,
    all_twelve_pressure_viscous_time_and_nonlinear_localization_terms_retained=True,
    both_background_pressure_tails_and_full_fixed_gauge_memory_retained=True,
    zero_initial_datum_joined_to_same_field=True,
    compact_smooth_global_force_relative_to_retained_lemmas_joined=True,
    one_Borel_sequence_controls_every_fixed_mixed_force_jet=True,
    whole_space_preterminal_energy_and_total_dissipation_bounds_joined=True,
    inherited_core_growth_preserved_on_localization_plateau=True,
    same_labels_cutoff_sequences_supports_width_and_parameters_preserved=True,
    exact_Gaussian_ledger_never_forward_inverted_or_newly_summed=True,
    historical_frozen_evidence_preserved_and_replayed=True,
    exact_polynomial_and_synthetic_norm_controls_are_not_continuum_proofs=True,
    retained_exact_core_and_finite_existence_inverse_gap_lemmas_still_dependencies=True,
    actual_infinite_coefficients_terminal_norms_or_Borel_scales_evaluated=False,
    analytic_or_zero_heat_extension_used=False,
    force_assumed_divergence_free=False,
    regular_velocity_or_pressure_at_singular_origin_inferred=False,
    future_global_velocity_or_pressure_constructed=False,
    global_competitor_pressure_recovery_and_comparison_joined=False,
    all_retained_analytic_dependencies_discharged=False,
    unconditional_admissible_force_complete_residual_or_blowup_verified=False,
    formal_certificate=False)
REMAINING=['recover_pressure_and_close_expanding_ball_comparison_for_every_admissible_global_competitor_of_this_same_zero_datum_and_compact_force_without_unproved_pressure_derivative_or_flux_assumptions',
    'discharge_retained_exact_core_and_finite_analytic_existence_inverse_gap_dependencies_before_unconditional_force_residual_or_blowup_claims']
DEPENDENCIES=('endpoint_foundation','realization_foundation','current_localized_force','localized_force_extension',
    'finite_induction_foundation','compact_mean_foundation','coefficient_support_audit','global_moment_stress')

def protocol():
    return dict(parent_commit='4899fa23a49f654c53bd3332c4c269bf9072e744',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        unchanged_parameters='Md64; unit viscosity; same h,N,positive lambda,width C_setup^-20,roots,geometry,labels,cutoff sequences and pressure gauge',
        domain='same current realization restricted to the existing late observation cylinder; cutoff support closure lies inside its domain',
        enlarged_observation='1/4<=q/Q<=4; q_e=min(qstar/16,2^(-ell0-2)); tau0=q_e/16; z0=(q_e/16)^D; 4Q<=16q<=6q_e<=3qstar/8',
        localization='c=chi_r(r^2)*chi_z(z)*chi_t(t); u=curl(c*A)+c*B=c*u0+grad(c) cross A; p=c*p0; full fixed gauge',
        derivative_orders='force C^m needs A C^(m+3), B C^(m+2), p0 C^(m+1); uniform next time jet reserves one additional order',
        pressure='-Delta p=partial_i partial_j(u_i*u_j)-div f; force need not be divergence free',
        corner='all fixed force jets are uniformly O(q^N) for every N; q<=2*tau+4*abs(z)^(1/D); all terminal force jets flat at origin',
        extension='sigma=t-1; f_plus=sum chi0(b_j*sigma)*sigma^j/j!*F_j; b0=2; increasing b_j=2^n_j; one sequence for a+m<=floor(j/2)',
        whole_support='||D_x^a partial_sigma^m T_j||<=C_jam*b_j^(m-j)<=2^(-j); tail<=2^(1-J), J>=max(1,2*(a+m))',
        force_support='K x [1-tau0,3/2]; smooth compact force on R3 x (0,infinity), zero initial datum',
        energy='E(t)+2*integral_0^t D <= (integral_0^t ||f||_2)^2; separately sup E<=F(1)^2 and total D<=F(1)^2/2; do not add separate extrema with same budget',
        growth='on inherited protected core and localization plateau: u_theta >= (e0/2)*tau^(-A), relative to retained core lemma',
        diagnostics='9 exact localization vector components at 3 positive viscosities; 12 nonzero omission controls; 98 Leibniz and 49 terminal factorial controls; exact compact polynomial energy integrals',
        scope='compact smooth global force and preterminal energy relative to retained exact-core and finite analytic lemmas; no global comparison or unconditional blow-up')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'endpoint_foundation';sys.path.insert(0,str(p))
    try:return load('force_foundation_current_endpoint_parent',p/'audit.py')
    finally:sys.path.pop(0)

@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'endpoint_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,d in f['source_hashes'].items():history['endpoint_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    return r,history

def sources():
    paths=[HERE/n for n in ('README.md','localized_bounds.py','localized_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'realization_foundation/realized_checks.py',HERE.parent/'current_localized_force/forcing.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():
    r=dict(current_axisymmetric_localization=domain_review(),global_uniform_terminal_force_jets=terminal_review(),
        one_whole_support_Borel_sequence=extension_review(),current_preterminal_energy_budget=energy_review(),admissibility_and_conditional_scope=scope_review())
    r.update(exact_reviews());return deepcopy(r)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews();gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(current_endpoint_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Current compact force scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Comparison or analytic dependencies erased')
    if set(r.get('gates',{}))!=set(reviews())|{'current_endpoint_parent_reproduced','historical_bytes_preserved'} or not all(r['gates'].values()):raise ValueError('Force foundation gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=997 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Compact smooth force and preterminal energy relative to retained lemmas; competitor comparison and unconditional claims remain open.')

if __name__=='__main__':main()
