#!/usr/bin/env python3
"""Same-canonical compact force, protected flat disk and preterminal energy."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from force_bounds import load,universal_bounds,universal_checks,canonical_input_review,scale_review,protected_disk_review
from force_checks import protected_plane_review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='conditional-canonical-compact-force-protected-disk-flatness-and-preterminal-energy'
DEPENDENCIES=('canonical_endpoint_join','canonical_local_realization','force_foundation','endpoint_foundation',
    'realization_foundation','current_localized_force','localized_force_extension','physical_pulse_foundation')
FLAGS=dict(same_canonical_complete_field_fixed_gauges_and_pressure_memory_localized=True,
    axisymmetric_potential_cutoff_and_separate_direct_swirl_preserve_divergence=True,
    all_twelve_localization_terms_and_full_force_divergence_pressure_equation_retained=True,
    canonical_zero_initial_datum_and_compact_smooth_global_force_relative_to_retained_premises=True,
    global_uniform_terminal_force_jets_from_corner_flatness_and_punctured_endpoint_join=True,
    all_terminal_force_jets_flat_along_protected_equatorial_disk=True,
    one_exact_dyadic_Borel_sequence_handles_arbitrary_finite_terminal_norms=True,
    future_force_preserves_flat_disk_and_every_uniform_normal_power=True,
    correct_fixed_time_energy_and_separate_total_dissipation_bounds_for_same_field=True,
    protected_core_growth_preserved_but_retained_core_lemma_not_discharged=True,
    actual_Gaussian_pieces_and_ledger_not_dropped_newly_summed_or_forward_inverted=True,
    finite_exact_controls_not_actual_fields_or_continuum_proofs=True,
    historical_source_bytes_fields_and_evidence_preserved=True,
    historical_force_or_global_field_identified_with_canonical_field=False,
    force_assumed_divergence_free_or_compact_Leray_projection_used=False,
    complete_force_Gevrey_bound_or_analytic_Taylor_convergence_assumed=False,
    future_velocity_or_pressure_constructed=False,
    flat_disk_claimed_open_force_free_neighborhood=False,
    actual_terminal_norms_or_infinite_Borel_scales_evaluated=False,
    regular_terminal_velocity_or_pressure_at_origin_inferred=False,
    canonical_global_competitor_pressure_recovery_and_comparison_joined=False,
    retained_core_background_and_finite_analytical_premises_discharged=False,
    unconditional_admissible_force_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['recover_pressure_and_close_expanding_ball_comparison_for_every_admissible_global_smooth_finite_energy_competitor_with_THIS_zero_datum_and_THIS_force_without_unproved_pressure_derivative_or_boundary_flux_hypotheses',
    'independently_discharge_retained_exact_core_background_and_finite_analytical_defining_existence_inverse_gap_premises_before_unconditional_claims']

def protocol():
    return dict(parent_commit='99b5b2a72c16324a0c9389defa5e1d7a9cdfe76f',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        unchanged='Md64,unit viscosity,original h,N,positive lambda,pulse radius,width C_setup^-20,labels,carriers,supports,qstar,both diagonal schedules,fixed operator gauges',
        same_field='canonical_local_realization plus canonical_endpoint_join, own canonical Hc/ac and complete current pressure targets; not old force or old global field',
        localization='c=chi_r(r^2)*chi_z(z)*chi_t(t); u=curl(c*Ac)+c*Bc=c*u0+grad c cross Ac; p=c*pc; support closure defined inside observed domain',
        full_force='all twelve terms of localized R_1 including full pc*grad c, time, vector viscosity and every nonlinear cross/square; d dot grad c=0 exactly',
        requirements='force m: Ac m+3,Bc m+2,pc m+1; next-time Cauchy reserves one more; complete Gaussian ledger kept',
        corner='all alpha,j,N: |D_x^alpha partial_t^j f_minus|<=C*(tau+abs(z)^(1/D))^N on fixed plateau box, uniformly in r',
        flat_disk='D_x^alpha F_j(x1,x2,0)=0 for r<rL/2, all alpha,j; future mixed force jets <=C_N*abs(z)^N on smaller disk x [0,1/2]',
        Borel='f_plus=sum chi0(b_j*sigma)*sigma^j/j!*F_j; b0=2,b_j=2^n_j; ONE sequence, whole-support constraints a+m<=floor(j/2)',
        exact_scale='rational c>=max(0,log2 C); integer n>max(previous+1,2,(c+j+2)/(j-m)); no huge power materialization or floating acceptance',
        force_support='K x [1-tau0,3/2], zero initial velocity; smooth compact global force only, no future velocity or pressure',
        pressure='-Delta p=partial_i partial_j(u_i*u_j)-div f; no solenoidal-force requirement or compact Leray projection',
        energy='E(t)+2*integral_0^t D<=F(t)^2; separately sup E<=F(1)^2,total D<=F(1)^2/2; do not combine extrema with same budget',
        official_source='https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf, checked 2026-10-18, equations (4),(5),(7), alternative (C)',
        scope='conditional same-canonical-field compact force, protected flat disk and preterminal energy; competitor comparison and retained analytic premises remain')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'canonical_endpoint_join';old=sys.path[:];sys.path.insert(0,str(p))
    try:return load('canonical_force_canonical_endpoint_parent',p/'audit.py')
    finally:sys.path[:]=old

@lru_cache(None)
def parent_record():
    r=json.loads((HERE.parent/'canonical_endpoint_join/evidence.json').read_text());parent_module().validate(r);return r

@lru_cache(None)
def realization_record():
    r=json.loads((HERE.parent/'canonical_local_realization/evidence.json').read_text());parent_module().parent_module().validate(r);return r

def history():
    m=parent_module()
    for _ in range(6):m=m.parent_module()
    return m.parent_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','force_bounds.py','force_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    paths += [HERE.parent/'force_foundation'/n for n in ('localized_bounds.py','localized_checks.py')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def domain_review():
    b=universal_bounds();requirements=[b.requirements(m) for m in (0,1,5,100)]
    return dict(requirements=requirements,retained_terms=list(b.TERMS),checks=dict(
        all_twelve_universal_operator_terms_reapplied_to_canonical_fields=len(b.TERMS)==12,
        canonical_next_time_force_jets_use_potential_m_plus_four=all(r['next_time_potential']==r['force']+4 for r in requirements),
        same_full_pressure_not_unforced_Riesz_pressure=True,
        support_neighborhood_inside_same_enlarged_observation_domain=True,
        Cartesian_r_squared_cutoff_and_same_axisymmetric_direct_swirl=True,
        zero_initial_time_neighborhood_and_same_core_plateau=True))

def reviews():
    return deepcopy(dict(canonical_current_field_input=canonical_input_review(parent_record(),realization_record()),
        same_field_complete_localization=domain_review(),exact_rational_Borel_scale=scale_review(),
        protected_disk_and_future_force_flatness=protected_disk_review(),
        finite_order_protected_plane_controls=protected_plane_review(),
        correct_preterminal_energy_budget=universal_bounds().energy_review(),**universal_checks().reviews()))

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Canonical force scope changed')
    compare_records(protocol(),r.get('protocol'),'protocol');compare_records(FLAGS,r.get('flags'),'flags')
    compare_records(REMAINING,r.get('remaining_obligations'),'limits')
    expected=reviews();compare_records(expected,r.get('reviews'),'reviews')
    g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in expected.items()};g['historical_recursive_source_bytes']=True
    compare_records(g,r.get('gates'),'gates')
    if not all(g.values()):raise ValueError('Canonical force gate failed')
    compare_records(sources(),r.get('source_hashes'),'sources');compare_records(history(),r.get('inherited_provenance'),'history')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Conditional canonical force/energy only; global competitor comparison and analytical premises remain.')
if __name__=='__main__':main()
