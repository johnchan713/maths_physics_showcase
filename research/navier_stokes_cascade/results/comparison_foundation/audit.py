#!/usr/bin/env python3
"""Current compact-force comparison with an improved R^-2 infinity remainder."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from compact_bounds import load,interface_review,pressure_review,radius_review,limit_review,scope_review
from compact_checks import reviews as exact_reviews

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-source-compact-pressure-comparison-relative-to-retained-core-and-finite-analytic-lemmas'
FLAGS=dict(same_current_localized_field_force_zero_datum_and_complete_pressure_used=True,
    current_force_and_energy_foundations_joined_to_comparison_hypotheses=True,
    pressure_gradient_recovered_without_assumed_competitor_pressure_decay=True,
    common_force_divergence_cancels_exactly_once_in_difference=True,
    compact_mixed_stress_pressure_estimated_outside_its_support=True,
    full_expanding_ball_energy_remainder_improved_from_four_thirds_to_two=True,
    no_global_competitor_derivative_or_dissipation_hypothesis_added=True,
    general_comparison_lemma_has_explicit_ordinary_hypotheses=True,
    current_global_smooth_finite_energy_competitor_exclusion_relative_to_retained_lemmas=True,
    same_parameters_width_labels_supports_gauge_and_diagonal_sequences_preserved=True,
    exact_Gaussian_ledger_never_forward_inverted_or_newly_summed=True,
    exact_local_polynomials_and_exterior_point_stresses_are_diagnostics=True,
    historical_frozen_evidence_preserved_and_replayed=True,
    retained_exact_core_and_finite_existence_inverse_gap_lemmas_still_dependencies=True,
    actual_infinite_field_or_global_competitor_evaluated=False,
    comparison_constants_uniform_through_t_one_assumed=False,
    pressure_Poisson_equation_alone_used_to_remove_harmonic_gradient=False,
    radius_minus_two_claimed_universal_optimal_rate_for_each_fixed_solution=False,
    future_velocity_or_pressure_constructed=False,
    all_retained_analytic_dependencies_discharged=False,
    unconditional_forced_breakdown_or_unforced_singularity_verified=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['trace_retained_exact_core_to_the_infinite_space_contraction_original_equations_and_identification_with_the_completed_global_backward_pressure',
    'independently_discharge_the_retained_finite_background_pulse_defining_existence_inverse_and_gap_lemmas_at_the_unchanged_parameters',
    'identify_and_close_the_first_unproved_continuum_identity_or_inequality_before_unconditional_theorem_assessment']
DEPENDENCIES=('force_foundation','global_solution_comparison','current_solution_comparison','endpoint_foundation',
    'realization_foundation','core_all_order_foundation','core_contraction_audit','finite_induction_foundation')

def protocol():
    return dict(parent_commit='44fc3767b90e07f11df62b4726d2d48a3d8e85ec',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        unchanged_parameters='Md64; unit viscosity; same h,N,positive lambda,width C_setup^-20,roots,geometry,labels,cutoff sequences and pressure gauge',
        current_input='force_foundation: same u=curl(c*A)+c*B, p=c*p0, zero datum, compact smooth global force and preterminal energy',
        competitor='same force/datum; smooth on compact sets; sup_[0,T] ||v||_2 finite for each T<1; no global pressure derivative or dissipation assumption',
        pressure='g=vv-uu=ww+h; h=w tensor u+u tensor w; -Delta(P-p)=DDg; pi_star=sum R_i R_j g_ij',
        recovery='L1 stress gives H^-2 pressure; time-averaged momentum gives H^-3 harmonic gradient, whose weighted L2 Fourier function supported at zero vanishes',
        cutoff='chi_R=phi(x/R)^8; E_R=int chi_R*|w|^2; D_R^2=int chi_R*|grad w|^2; B_R=||phi_R^4*w||_6',
        compact_cross='support(h) subset K subset B_L; R>=max(1,2L); |pi_h|<=16/(4*pi)*R^-3*sum ||h_ij||_1 on cutoff annulus; flux O(R^-5/2)',
        quadratic_pressure='main L^(3/2) norm <=C*B_R*W_T; commutator L^(4/3) norm <=C*R^-3/4*W_T^2; pressure flux <=C_T*(R^-1*B_R^3/2+R^-7/4*B_R^3/4)',
        energy='0.5*E_R_prime+0.5*D_R^2<=||grad u||_infinity*E_R+C_T*R^-2; E_R(0)=0',
        limits='fix T<1, then R->infinity; E_R<=2*C_T*T*exp(2*T*L_T)*R^-2; only afterward approach the protected terminal core',
        historical='older valid R^-4/3 estimate preserved; no historical localized field substituted for the current one',
        external_input='Riesz Lp boundedness, Knapp Advanced Real Analysis III.5 Theorem 3.26 and III.6 Riesz example, printed pp83,93; checked 2026-10-17',
        source='https://www.math.stonybrook.edu/~aknapp/download/a2-1-realanal-clickable.pdf',
        diagnostics='9 exact stress identities, 3 complete local weighted-energy identities, 60 compact exterior tensor-kernel controls and exact radius exponents; not actual whole-space flows',
        scope='same-current-field global competitor exclusion relative to retained core and finite analytic lemmas; no unconditional breakdown theorem')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'force_foundation';sys.path.insert(0,str(p))
    try:return load('comparison_foundation_current_force_parent',p/'audit.py')
    finally:sys.path.pop(0)

@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'force_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,d in f['source_hashes'].items():history['force_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    return r,history

@lru_cache(None)
def general_review():
    p=HERE.parent/'global_solution_comparison';sys.path.insert(0,str(p))
    try:m=load('comparison_foundation_general_pressure_lemma',p/'audit.py')
    finally:sys.path.pop(0)
    r=m.generate();m.validate(r);f=json.loads((p/'evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    return dict(ordinary_general_lemma_not_historical_special_field=True,
        checks=dict(all_general_controls_reproduced=all(r['gates'].values()),
            general_source_bindings_unchanged=r['source_hashes']==f['source_hashes'],
            all_general_inherited_sources_preserved=all(r['inherited_provenance'].values()),
            current_refinement_proved_separately_not_claimed_in_old_record=True))

def sources():
    paths=[HERE/n for n in ('README.md','compact_bounds.py','compact_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'realization_foundation/realized_checks.py',HERE.parent/'global_solution_comparison/comparison_bounds.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():
    r=dict(current_force_and_comparison_interface=interface_review(parent_review()[0]),
        pressure_gradient_and_compact_mixed_stress=pressure_review(),improved_complete_radius_budget=radius_review(),
        fixed_interval_then_protected_core_limit=limit_review(),conditional_scope_and_foundation_gate=scope_review(),
        general_pressure_comparison_controls=general_review())
    r.update(exact_reviews());return deepcopy(r)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews();gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(current_force_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Current comparison scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Foundation obligations erased')
    if set(r.get('gates',{}))!=set(reviews())|{'current_force_parent_reproduced','historical_bytes_preserved'} or not all(r['gates'].values()):raise ValueError('Comparison foundation gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=1015 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Radius remainder R^-2 and same-current-field comparison relative to retained analytic lemmas; unconditional breakdown remains unverified.')

if __name__=='__main__':main()
