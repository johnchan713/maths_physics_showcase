#!/usr/bin/env python3
"""Exact reference angular/heat roots and the scalar terminal waiting root."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from root_bounds import load,matrix_review,target_review,smallness_review,waiting_review,scope_review
from root_checks import reviews as exact_reviews

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='reference-angular-heat-and-wait-roots-audited-relative-to-remaining-construction-lemmas'
FLAGS=dict(actual_reference_angular_and_heat_moment_maps_traced=True,
    uniform_exact_bump_matrix_bounds_proved_without_quadrature=True,
    full_C1_value_plus_derivative_norm_retained=True,
    sharper_angular_inverse_three_and_heat_inverse_fifteen_supplied=True,
    angular_root_C1_bound_4200_lambda30_supplied=True,
    heat_root_C1_bound_1920_h4_over_Xstar_supplied=True,
    same_original_accepted_roots_identified_by_uniqueness=True,
    scalar_terminal_wait_positive_finite_unique_and_eta_independent=True,
    all_fixed_one_sided_eta_jets_of_same_roots_supplied=True,
    exact_pressure_and_energy_and_renormalized_I_rows_retained=True,
    finite_h_infinite_angular_integral_not_erased=True,
    same_parameters_width_supports_labels_gauge_and_diagonal_schedules_preserved=True,
    historical_frozen_evidence_preserved_and_replayed=True,
    rational_bounds_are_continuum_enclosures_finite_quadratures_are_diagnostics=True,
    actual_selected_lambda_coefficients_or_full_profile_evaluated=False,
    matrix_parameter_range_promoted_to_target_range=False,
    constant_Q_rule_applied_to_eta_dependent_annular_map=False,
    common_complex_strip_for_heat_targets_claimed=False,
    axial_two_moment_and_scalar_strength_dependencies_discharged=False,
    all_annular_compact_repair_and_pulse_existence_inverse_gap_dependencies_discharged=False,
    complete_core_global_pressure_identification_unconditionally_verified=False,
    all_retained_analytic_dependencies_discharged=False,
    unconditional_force_residual_or_blowup_verified=False,formal_certificate=False)
REMAINING=['independently_verify_the_actual_axial_two_moment_solve_and_scalar_strength_root_including_all_later_energy_and_nonzero_end_bump_terms',
    'verify_annular_and_compact_repair_targets_inverses_and_gaps_from_actual_maps_and_complete_the_core_global_backward_pressure_identification',
    'discharge_remaining_finite_background_and_pulse_defining_existence_inverse_gap_lemmas_before_unconditional_theorem_assessment']
DEPENDENCIES=('comparison_foundation','pulse_moment_audit','heat_exterior_audit','global_moment_stress',
    'annular_exterior_foundation','core_all_order_foundation','modulation_repair_all_order','axis_matching_audit')

def protocol():
    return dict(parent_commit='bf4794ad6e0dc309acab0cc6e84ad74a8ea558a1',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        unchanged_parameters='Md64; unit viscosity; same h,N,positive lambda,width C_setup^-20,roots,geometry,labels,cutoff sequences and pressure gauge',
        actual_parameters='T=exp(64)+10>=64; lambda=exp(-4T)<=exp(-256); h=lambda^2; Xstar>=100; co=1/1000',
        bump='beta_j(y)=sigma_prime((y-center_j)/(3/10)+1/2); exact mass 3/10; sup sigma_prime<=64; sup abs sigma_second<=20000',
        angular='centers 0,2 in local first-bump coordinates; rates 1-lambda,-1-2lambda; pressure row factor two and quadratic row factor one',
        heat='centers 1/2,5/2,9/2; rates -1-2lambda,-2lambda,1-lambda; pressure and minus-S quadratic factor 1/2; I row linear',
        uniform_matrix='0<=lambda<=1/100; positive masses and Lagrange polynomial inverse; inverse norms <3 angular, <15 heat; quadratic norms <32 and <29',
        targets='angular C1<=700*lambda^30 only on 0<lambda<=exp(-256); heat C1<=64*h^4/Xstar from all three complete infinite heat integrals',
        Banach_norm='max_j(sup abs c_j+sup abs partial_eta c_j); constant B and Q; contraction requires 8*beta^2*q*delta<1',
        roots='angular max C1<=4200*lambda^30; heat max C1<=1920*h^4/Xstar; subsets of original uniqueness balls',
        jets='same roots; pointwise Jacobian inverse <=6 angular, <=30 heat; d_k<=2beta*(z_k+q*sum_(j=1)^(k-1) d_j*d_(k-j))',
        waiting='q_w=(1-h)*I/(XH)-1; omega=log(q_before/qp)/(1-h); qp=rho/(1-rho)*integral exp((1-h)*(1+2v))*sigma_prime(v) dv; rho=h/1000',
        pressure='angular Cp change zero; heat Cp, S and renormalized I changes cancel heat edit exactly as eta functions; other exact roots still required for full Pi_axis+Cp_infinity=0',
        diagnostics='44 exact inverse-product components, exact disjoint relative-square moment controls, 12 Catalan jet orders, terminal integrating-factor control and six moderate-parameter bump-integral comparisons',
        scope='two reference correction roots and scalar waiting root; remaining axial strength, annular/repair and finite pulse lemmas explicit; no unconditional construction')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'comparison_foundation';sys.path.insert(0,str(p))
    try:return load('root_foundation_current_comparison_parent',p/'audit.py')
    finally:sys.path.pop(0)

@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'comparison_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,d in f['source_hashes'].items():history['comparison_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    return r,history

def sources():
    paths=[HERE/n for n in ('README.md','root_bounds.py','root_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'realization_foundation/realized_checks.py',HERE.parent/'outer_pressure_pilot/schedule.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():
    r=dict(actual_uniform_bump_matrix_and_quadratic_bounds=matrix_review(),complete_angular_and_heat_targets=target_review(),
        same_root_smallness_and_strict_gaps=smallness_review(),exact_positive_scalar_wait=waiting_review(),remaining_foundation_scope=scope_review())
    r.update(exact_reviews());return deepcopy(r)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews();gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(current_comparison_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Reference root scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Other root or construction obligations erased')
    if set(r.get('gates',{}))!=set(reviews())|{'current_comparison_parent_reproduced','historical_bytes_preserved'} or not all(r['gates'].values()):raise ValueError('Reference root gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=1033 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Angular/heat exact roots and scalar wait supplied; other construction dependencies and unconditional breakdown remain open.')

if __name__=='__main__':main()
