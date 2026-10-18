#!/usr/bin/env python3
"""Canonical real covariance and fresh initial means, with strict scope guards."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from real_covariance import complete_review,variance_review
from mean_checks import pressure_review,order_review,nonflat_review,recompute_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='conditional-canonical-complete-real-covariance-and-fresh-initial-mean-join'
DEPENDENCIES=('canonical_primary_foundation','compact_mean_foundation','complete_covariance_initial_means',
              'physical_wave_interface','harmonic_closure_foundation')
FLAGS=dict(canonical_real_complete_curl_all_six_covariance_inputs_joined=True,
    real_primary_sine_cosine_cancellation_gives_quadratic_curl_square=True,
    canonical_Haar_principal_target_matches_with_its_own_rebuilt_roots=True,
    pressure_fast_recompute_slow_recompute_prefix_uses_current_canonical_wave=True,
    conditional_initial_three_defects_S2_and_tangential_mean_M2_minus_3kappa=True,
    axial_velocity_full_pressure_moment_moving_bump_and_both_base_tails_retained=True,
    exact_counterexample_two_fluxes_do_not_fix_tangential_variance=True,
    polynomial_vs_epsilon_obstruction_prevents_automatic_mean_target_transfer=True,
    every_fixed_jet_uses_same_value_domain_and_retained_smooth_support_hypotheses=True,
    historical_fields_and_frozen_evidence_preserved=True,
    exact_finite_controls_distinguished_from_continuum_estimates=True,
    two_selected_fluxes_determine_all_six_covariances=False,
    old_mean_targets_or_old_mean_velocities_transplanted=False,
    inverse_length_closeness_is_extra_epsilon_power=False,
    real_primary_cancellation_applies_to_generic_complex_later_wave=False,
    nonlinear_defects_identically_zero_after_five_linear_rows=False,
    actual_canonical_family_saturates_counterexample_lower_bound=False,
    all_later_signed_sourced_maps_and_complete_residual_closed=False,
    unconditional_pulse_family_or_blowup_proved=False,formal_certificate=False)
REMAINING=['apply_canonical_signed_covariance_inverse_and_recompute_complete_complex_wave_means_using_canonical_roots',
    'join_every_generated_supported_source_and_finite_cycle_on_same_current_canonical_family',
    'full_physical_scale_viscous_cutoff_interaction_residual_with_axial_velocity_and_complete_pressure_memory',
    'retained_exact_core_background_and_global_construction_premises_before_unconditional_claims']


def protocol():
    return dict(parent_commit='d3c4dd1b00337d54d446888355802a9d8a824bf8',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        construction='canonical reference-mixing real primary with its own rebuilt positive roots; same Md64,h,N,lambda>0,carriers,r0,width C^-20,gauge,unit viscosity',
        full_curl='w=sqrt(epsilon)*(A cos(kPhi)-B sin(kPhi)/k); B=curl_*(n cross A/|n|^2); all cutoff and cylindrical derivatives',
        real_covariance='W=epsilon/2*(A tensor A+B tensor B/k^2); all six entries; mixed quadratures vanish only for real primary',
        canonical_inputs='W in M1; selected Haar flux error M_(2-2kappa); full real curl square not generic linear cross',
        obstruction='same transverse plane and rebuilt positive weights can leave eta/2*(0,1,-1) tensor (0,1,-1); epsilon/L is not epsilon^(1+delta)*S^B',
        order='pressure from zero means; fast zero-Haar velocity; recompute full pressure and three defects; slow fixed five-map on current targets; recompute again',
        initial_result='three defects S2, radial mean M2 and full tangential angular-mean residual M_(2-3kappa), plus retained flat terms; not nonzero-angular residual',
        moment='Jtheta=int R^2*(Gv+Vgamma+gamma*v+Wtheta z); Jz=int R*(2Ggamma+gamma^2+Wzz)-int R^2*gr/2',
        memory='int R*bar(pm)=-int R^2*bar(gr)/2+c_rho*P; whole partial_Z(Jz+c_rho*P); current axial velocity; same complete Pi1/Pi2 base subtraction',
        support='same accepted enlarged charts and full I4 potential hull; actual A0/A1 velocities retained in products; Gaussian cutoff tails additive not inverse sources',
        diagnostic_scope='exact Fraction quadrature and polynomial means; finite endpoint regularity and plane-wave counterexample not actual giant canonical realization',
        remaining='later canonical signed/sourced complex states and full residual require fresh consistent construction; no unconditional theorem')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'canonical_primary_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('canonical_mean_current_primary',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old


def history():
    harmonic=parent_module().parent_module().parent_module().parent_tools()
    return harmonic.verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])


def sources():
    files=[HERE/n for n in ('README.md','audit.py','real_covariance.py','mean_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def reviews():return deepcopy(dict(canonical_complete_real_curl=complete_review(),unfixed_tangential_variance=variance_review(),
    fresh_pressure_moments_and_targets=pressure_review(),canonical_initial_order_ledger=order_review(),
    polynomial_not_epsilon_gain=nonflat_review(),initial_Haar_and_nonlinear_recompute=recompute_review()))


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
    if r.get('status')!=STATUS:raise ValueError('Canonical initial scope changed')
    compare_records(protocol(),r.get('protocol'),'protocol');compare_records(FLAGS,r.get('flags'),'flags')
    compare_records(REMAINING,r.get('remaining_obligations'),'limits')
    expected=reviews();compare_records(expected,r.get('reviews'),'reviews')
    g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in expected.items()};g['historical_recursive_source_bytes']=True
    compare_records(g,r.get('gates'),'gates')
    if not all(g.values()):raise ValueError('Canonical initial gate failed')
    compare_records(sources(),r.get('source_hashes'),'sources');compare_records(history(),r.get('inherited_provenance'),'history')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if a.verify_record:
        old=json.loads(a.verify_record.read_text());validate(old);compare_records(r,old)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Conditional canonical INITIAL mean join only; later complex cycles and full residual remain.')
if __name__=='__main__':main()
