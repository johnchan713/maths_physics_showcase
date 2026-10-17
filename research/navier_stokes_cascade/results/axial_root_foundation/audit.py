#!/usr/bin/env python3
"""Actual reference axial moment and complete scalar-strength foundation."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from axial_bounds import HERE,load,matrix_review,shape_review,budget_review,prefix_tail_review,scope_review
from axial_checks import reviews as exact_reviews

PROJECT=HERE.parent.parent
STATUS='reference-axial-moment-and-complete-strength-roots-audited-relative-to-remaining-construction-lemmas'
FLAGS=dict(actual_axial_two_moment_map_and_complete_strength_equation_traced=True,
    uniform_scaled_inverse_below_four_supplied=True,
    full_fixed_trial_C1_coefficients_below_lambda20=True,
    actual_composed_C1_coefficients_below_two_lambda20=True,
    complete_prefix_and_all_post_energy_stages_bounded=True,
    nonzero_end_bump_quadratic_linear_and_constant_terms_retained=True,
    positive_unique_strength_in_same_frozen_bracket=True,
    root_Jacobian_above_point_four_nine_supplied=True,
    root_eta_derivative_below_Tlambda_over_thirty_supplied=True,
    all_fixed_one_sided_root_and_coefficient_jets_supplied=True,
    heat_and_compensation_exact_net_S_identity_retained=True,
    reference_power_identity_distinguished_from_final_heat_tail=True,
    historical_frozen_evidence_preserved_and_replayed=True,
    same_parameters_roots_supports_width_labels_gauge_and_cutoff_sequences=True,
    actual_selected_lambda_or_huge_support_coordinates_evaluated=False,
    matrix_invertibility_at_lambda_zero_claimed=False,
    strength_evenness_assumed=False,
    finite_quadrature_diagnostics_promoted_to_continuum_proof=False,
    annular_and_compact_repair_dependencies_discharged=False,
    complete_core_global_pressure_identification_unconditionally_verified=False,
    all_retained_analytic_dependencies_discharged=False,
    unconditional_complete_residual_or_blowup_verified=False,formal_certificate=False)
REMAINING=['audit_actual_annular_and_compact_repair_targets_inverses_gaps_and_eta_dependent_bilinear_terms',
    'complete_the_current_core_global_backward_pressure_identification_with_all_matching_roots_and_gauge_retained',
    'discharge_remaining_finite_background_and_pulse_existence_inverse_gap_lemmas_before_unconditional_theorem_assessment']
DEPENDENCIES=('root_foundation','pulse_moment_audit','heat_exterior_audit','global_moment_stress',
    'annular_exterior_foundation','core_all_order_foundation','modulation_repair_all_order','axis_matching_audit')

def protocol():
    return dict(parent_commit='a41b941add8dc9c0368e3cea7c4f0aed979a1f78',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        parameters='same Md64,T=exp(64)+10,lambda=exp(-4T),h=lambda^2,Tw=240T,Tf=1000,co=1/1000,width C_setup^-20 and unit viscosity',
        matrix='B_ij=b(si)*exp(2si*j); s1=1/2-lambda,s2=1/2-2lambda; positive bump mass 3/10, width 3/10',
        inverse='0<lambda<=1/100; scaled infinity inverse norm lambda*||B^-1||<4; lambda=0 singular',
        targets='0<=Amp<=6/5; full eta C1 norm and strength derivative <=200exp(-4/(5lambda)); coefficient bound 800exp(-4/(5lambda))/lambda<=lambda^20',
        prefix='F1<=11/6 independently; Ka,Ja<=4,K2a<=16; full normalized s0<=500T and abs s0_eta<1',
        continuation='all eight later stages and signed angular edits; R<=256T and abs R_eta<=1024T; terminal waiting and nonzero end bump retained',
        scalar='F=Amp^2*K-(1-exp(-26))/4+lambda*s0-lambda*exp(-26)*R/2+lambda*sum(wj*cj(Amp,eta)^2)',
        shape='1024 positive monotone boxes, 80 interval digits; 0.2450496182718<K<0.2450496221363; no quadrature-based enclosure',
        roots='same 1.0100502<Amp<1.0100504; unique positive root, Jacobian>49/100, abs Amp_eta<Tlambda/30; all fixed one-sided jets',
        composition='actual coefficients max component full eta C1 norm<=2lambda^20; fixed trial bound is distinct',
        heat='complete net pressure, energy and renormalized-I correction identities retained; S-strength equation unchanged by the heat attachment and compensation',
        diagnostics='exact polynomial two-moment and affine energy identities; moderate-parameter actual bump matrices and stable solver; proxy scalar discriminant controls',
        scope='reference axial two-moment and scalar strength roots; annular/repair, global pressure identification and finite-pulse lemmas still required')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'root_foundation';sys.path.insert(0,str(p))
    try:return load('axial_root_current_root_parent',p/'audit.py')
    finally:sys.path.pop(0)

@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'root_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,d in f['source_hashes'].items():history['root_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    return r,history

def sources():
    paths=[HERE/n for n in ('README.md','axial_bounds.py','axial_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'pulse_moment_audit/pulse.py',HERE.parent/'pulse_moment_audit/bounds.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():
    r=dict(actual_scaled_axial_inverse=matrix_review(),positive_shape_enclosure=shape_review(),
        complete_scalar_root_and_derivative_budget=budget_review(),independent_prefix_and_complete_tail=prefix_tail_review(),
        remaining_foundation_scope=scope_review());r.update(exact_reviews());return deepcopy(r)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews();gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(current_root_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Axial root scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Other construction obligations erased')
    if set(r.get('gates',{}))!=set(reviews())|{'current_root_parent_reproduced','historical_bytes_preserved'} or not all(r['gates'].values()):raise ValueError('Axial root gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=1051 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Reference axial and strength roots supplied; remaining matching maps and unconditional construction remain open.')

if __name__=='__main__':main()
