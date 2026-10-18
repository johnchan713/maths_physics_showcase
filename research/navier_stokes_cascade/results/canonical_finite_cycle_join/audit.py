#!/usr/bin/env python3
"""Fresh canonical finite-cycle bridge; no historical field identification."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from cycle_bounds import input_review,prefix_review,signed_review,cycle_review,physical_review,diagonal_review
from freshness import review as freshness_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='conditional-canonical-seed-and-all-finite-cycle-operator-bridge'
DEPENDENCIES=('canonical_initial_mean_join','canonical_primary_foundation','regular_source_foundation',
    'first_sourced_foundation','current_complex_foundation','finite_induction_foundation',
    'sourced_inverse_foundation','harmonic_closure_foundation')
FLAGS=dict(fresh_canonical_initial_state_rebuilt_through_signed_sourced_and_mean_prefix=True,
    one_fixed_canonical_H_and_positive_primary_roots_for_all_signed_maps=True,
    conditional_canonical_seed_B0_three_halves_minus_two_kappa_supplied=True,
    universal_operator_estimates_reapplied_to_current_canonical_inputs=True,
    conditional_all_finite_cycle_gain_half_minus_kappa_on_whole_half_line=True,
    current_supported_source_expressions_and_complete_residual_groups_exhausted=True,
    complete_pressure_viscosity_all_six_covariances_and_nonlinear_squares_retained=True,
    recomputation_after_every_wave_and_mean_change_and_current_postfast_targets=True,
    same_full_base_Pi1_Pi2_axial_velocity_and_whole_moving_pressure_moment=True,
    old_additive_ledger_preserved_new_tails_appended_never_forward_inverted=True,
    fixed_finite_physical_target_and_jet_has_a_conditional_finite_state=True,
    all_fixed_jets_same_value_domain_and_original_supports_width_and_carriers=True,
    exact_controls_and_metadata_model_not_continuum_or_actual_field_certificates=True,
    historical_pulses_sources_roots_mean_fields_or_realization_identified_as_canonical=False,
    old_source_or_mean_targets_reused_by_polynomial_closeness=False,
    initial_real_covariance_shortcut_used_for_complex_cycle=False,
    stage_constants_uniform_or_fixed_q_norm_contraction_established=False,
    arbitrary_stage_j_of_q_substitution_licensed=False,
    one_simultaneous_canonical_infinite_realization_constructed=False,
    retained_core_background_and_all_analytic_premises_discharged=False,
    unconditional_complete_residual_improvement_or_blowup_proved=False,formal_certificate=False)
REMAINING=['construct_one_consistent_canonical_complete_potential_realization_with_controlled_diagonal_cutoffs_direct_swirl_and_pressure_gauges',
    'join_its_endpoint_localization_force_energy_and_global_comparison_to_that_same_field',
    'discharge_retained_exact_core_background_and_finite_analytic_gap_premises_before_unconditional_claims']


def protocol():
    return dict(parent_commit='5e49006bab195749c75635a710e1d13216ee7f2b',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        family='canonical reference-mixing primary, own Hc and positive ac fixed for every stage; historical families not identified',
        fixed='Md64,unit viscosity,original h,N,lambda>0,r0,width C_setup^-20,carriers,common qstar and closed enlarged supports',
        prefix='fresh real signed; current fast/slow means; first regular sourced; current generic signed; fresh current fast/slow means; recompute after every change',
        seed='B0=3/2-2kappa; full tangential mean B0+1/2; defects 5/2-kappa; rho=mu=1-kappa; accumulated Haar tangential 1 radial 2',
        cycle='current regular sourced wave; recompute; canonical generic signed; recompute; zero-Haar fast; recompute; current slow five-map; recompute',
        finite='Bj=3/2-2kappa+j*(1/2-kappa); Cj=Bj+1/2; all finite j, not uniform stage constants',
        gain_limit='gain<=min(1/2-kappa,rho-1/2,mu-1/2); stale postfast targets leave gap 2kappa',
        inverse='pi_m=i(n.K*t_m-nprime.t_m+n.f_m)/(km|n|^2); current sources, zero entrance; same operator not same old source',
        pressure='all six covariances,current axial velocities,Pi1/Pi2 and whole partial_Z(Jz+c_rho*P); actual A0/A1 velocities in products',
        tails='parent additive ledger retained; each new (1-psi)f+Tpsi appended exactly; never forward-inverted',
        physical='fixed j,m only: |R(Uc_j,pc_j)|m<=Cjm*q^(hBj-(2A+3/2+2m)); Mj<=4*2^j',
        diagonal_limit='2^(j^2)*q^(hBj-K0) satisfies every fixed-stage estimate but can grow for j=ceil(2h*gain*ell)',
        controls='exact Fraction signed maps, universal affine verifier, independent complete-curl mode fixtures, source-stamp model; no giant actual fields',
        scope='conditional finite-cycle consistency bridge, not unconditional pulse/residual theorem or one limiting realization')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'canonical_initial_mean_join';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('canonical_cycle_current_initial',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old


@lru_cache(None)
def initial_record():
    r=json.loads((HERE.parent/'canonical_initial_mean_join/evidence.json').read_text())
    parent_module().validate(r);return r


def history():
    m=parent_module()
    for _ in range(3):m=m.parent_module()
    return m.parent_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])


@lru_cache(None)
def complete_controls():
    p=HERE.parent/'finite_induction_foundation/finite_checks.py'
    spec=importlib.util.spec_from_file_location('canonical_cycle_exact_modes',p)
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m


def sources():
    files=[HERE/n for n in ('README.md','audit.py','cycle_bounds.py','freshness.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def reviews():
    return deepcopy(dict(canonical_current_seed_inputs=input_review(initial_record()),fresh_canonical_prefix=prefix_review(),
        fixed_canonical_signed_map=signed_review(),whole_half_line_operator_bridge=cycle_review(),
        fixed_finite_physical_orders=physical_review(),arbitrary_diagonal_counterexample=diagonal_review(),
        canonical_current_source_and_target_protocol=freshness_review(),**complete_controls().reviews()))


def compare_records(a,b,path='root'):
    if type(a) is not type(b):raise ValueError('Record type changed at '+path)
    if isinstance(a,dict):
        if a.keys()!=b.keys():raise ValueError('Record keys changed at '+path)
        for k in a:compare_records(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if len(a)!=len(b):raise ValueError('Record length changed at '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare_records(x,y,path+'.'+str(i))
    elif a!=b:raise ValueError('Exact record changed at '+path)


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Canonical finite scope changed')
    compare_records(protocol(),r.get('protocol'),'protocol');compare_records(FLAGS,r.get('flags'),'flags')
    compare_records(REMAINING,r.get('remaining_obligations'),'limits')
    expected=reviews();compare_records(expected,r.get('reviews'),'reviews')
    g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in expected.items()};g['historical_recursive_source_bytes']=True
    compare_records(g,r.get('gates'),'gates')
    if not all(g.values()):raise ValueError('Canonical cycle gate failed')
    compare_records(sources(),r.get('source_hashes'),'sources');compare_records(history(),r.get('inherited_provenance'),'history')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if a.verify_record:
        old=json.loads(a.verify_record.read_text());validate(old);compare_records(r,old)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Conditional fixed-finite canonical bridge; one realization and retained analytical premises remain.')
if __name__=='__main__':main()
