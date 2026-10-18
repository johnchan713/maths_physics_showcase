#!/usr/bin/env python3
"""Explicit Gaussian/Haar/tail constants, with guarded conditional scope."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from gaussian_bounds import curvature_review,barrier_review,tail_review,degeneracy_review
from gaussian_checks import integral_review,ode_mass_review,exact_haar_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-fixed-tilt-gaussian-haar-and-tail-constants-from-spectral-floors'
DEPENDENCIES=('harmonic_closure_foundation','leading_frame_foundation','primary_covariance_foundation',
              'physical_pulse_foundation','geometric_pulse_replay')
FLAGS=dict(explicit_whole_pulse_Gaussian_curvature_from_current_spectral_floors=True,
    actual_Riccati_radial_amplitude_lower_and_upper_factors_quantified=True,
    primary_Haar_mass_and_direction_concentration_keep_radius_and_angular_half=True,
    explicit_cutoff_tail_a0_equals_c_r0_over_25=True,
    all_fixed_tail_jets_flat_on_same_accepted_band_domain=True,
    fixed_positive_tilt_required_and_zero_tilt_limit_not_uniform=True,
    full_moving_normal_source_pressure_and_background_memory_retained=True,
    original_value_acceptance_inequalities_not_new_derivative_dependent_bands=True,
    exact_rational_and_manufactured_numerical_controls_separate_from_analysis=True,
    every_actual_frame_column_constant_numerically_materialized=False,
    actual_giant_first_band_or_profile_numerically_evaluated=False,
    positive_pointwise_lower_bound_for_complete_oscillating_velocity=False,
    uniform_Gaussian_curvature_as_tilt_tends_to_zero=False,
    all_core_and_complete_residual_dependencies_discharged=False,
    unconditional_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['current_actual_frame_and_column_comparison_constants_in_original_value_acceptance_budget',
    'dependency_explicit_core_roots_inverse_complete_covariance_fresh_mean_pressure_and_full_viscous_cutoff_interaction_residual_join',
    'retained_exact_core_and_complete_construction_dependencies_before_any_unconditional_proof']


def protocol():
    return dict(parent_commit='a70298fb0db13af479b2547b52e5bf1454b99e6d',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        fixed_choices='Md64; original h,N,lambda>0,positive tilt,roots,carriers,width C^-20,r0=1/6057072,supports,gauge,unit viscosity',
        spectral='lambda_star/16 <= lambda0 <= 16*lambda_upper on same enlarged charts',
        curvature='c=lambda_minus*u^2/2*[(1+9u^2/4)^-2+2*(1+u^2)^-2]; C=3*lambda_plus*u^2/2*[(1+u^2/4)^-1+2*(1+u^2)^-1]',
        amplitude='D=Cd+9CE/8; xminus=(7/8)exp(-D); xplus=(9/8)exp(D); original L>=8CE/mu',
        Haar='h_sigma=r0^2/sqrt(L)*j_sigma; jminus=Jmin*exp(-2/3)*xminus^2*exp(-C); jplus=Jmax*xplus^2*sqrt(pi/c)',
        concentration='Cdir=D0+D1*xplus^2*exp(C)/(2*xminus^2*c); L>=25 and original cone acceptance length',
        tail='P<=exp(-cL/50)<=exp(-a0*ell^2); a0=c*r0/25; square completion controls every fixed jet and target on same domain',
        support='original cutoff plateau 1/5 and support inside 1/3, safely contained in older 2/5 enlargement; no support change',
        memory='nprime and normal-source pressure; Pi1/Pi2; current axial velocity and whole moving pressure moment; exact additive Gaussian tails never forwarded',
        scope='conditional Gaussian/Haar/tail input estimate, not all actual frame constants or complete residual, infinite summability or blow-up')


@lru_cache(None)
def parent_tools():
    path=HERE.parent/'harmonic_closure_foundation';before=sys.path[:];sys.path.insert(0,str(path))
    try:
        spec=importlib.util.spec_from_file_location('gaussian_current_harmonic_provenance',path/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=before


def history():return parent_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])


def sources():
    files=[HERE/n for n in ('README.md','audit.py','gaussian_bounds.py','gaussian_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def reviews():return deepcopy(dict(whole_pulse_curvature=curvature_review(),actual_amplitude_barrier=barrier_review(),
    explicit_all_band_tail_rate=tail_review(),fixed_tilt_and_underflow_limits=degeneracy_review(),
    independent_closed_integral=integral_review(),noncommuting_frame_and_actual_cutoff_mass=ode_mass_review(),
    exact_Haar_area_and_angular_normalization=exact_haar_review()))


def compare_records(a,b,path='root'):
    if type(a) is not type(b):raise ValueError('Record type changed at '+path)
    if isinstance(a,dict):
        if a.keys()!=b.keys():raise ValueError('Record keys changed at '+path)
        for k in a:compare_records(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if len(a)!=len(b):raise ValueError('Record length changed at '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare_records(x,y,path+'.'+str(i))
    elif isinstance(a,float):
        if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-10):raise ValueError('Numerical diagnostic changed at '+path)
    elif a!=b:raise ValueError('Record changed at '+path)


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['recursive_frozen_bytes_preserved']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Gaussian scope changed')
    compare_records(protocol(),r.get('protocol'),'protocol');compare_records(FLAGS,r.get('flags'),'flags')
    compare_records(REMAINING,r.get('remaining_obligations'),'limits')
    expected=reviews();compare_records(expected,r.get('reviews'),'reviews')
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in expected.items()};gates['recursive_frozen_bytes_preserved']=True
    compare_records(gates,r.get('gates'),'gates')
    if not all(gates.values()):raise ValueError('Gaussian estimate gate failed')
    compare_records(sources(),r.get('source_hashes'),'sources');compare_records(history(),r.get('inherited_provenance'),'history')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if a.verify_record:
        old=json.loads(a.verify_record.read_text());validate(old);compare_records(r,old)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Conditional analytical Gaussian/Haar/tail bounds; manufactured controls are not a continuum proof.')
if __name__=='__main__':main()
