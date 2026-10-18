#!/usr/bin/env python3
"""Guarded canonical-primary reconstruction and quantitative root stability."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from canonical_bounds import budget_review,matrix_review,acceptance_review
from canonical_checks import physical_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='conditional-canonical-primary-rebuild-and-relative-root-stability'
DEPENDENCIES=('frame_column_foundation','gaussian_mass_foundation','primary_covariance_foundation',
              'physical_pulse_foundation','geometric_pulse_replay')
FLAGS=dict(canonical_manuscript_reference_frame_and_growing_entrance_used=True,
    canonical_Riccati_mass_direction_and_positive_principal_roots_derived_conditionally=True,
    same_operator_distinct_entrances_have_Gaussian_relative_difference_bound=True,
    full_normal_pressure_difference_and_original_Haar_area_retained=True,
    explicit_relative_root_stability_requires_common_value_acceptance=True,
    canonical_principal_roots_rebuilt_not_reused_from_inherited_frame=True,
    old_pulse_choices_and_frozen_evidence_preserved=True,
    fixed_jets_use_same_positive_value_domain_and_retained_flat_edge_hypotheses=True,
    exact_rational_and_numerical_diagnostics_separate_from_continuum_argument=True,
    canonical_and_inherited_physical_primaries_or_roots_exactly_equal=False,
    small_diagnostic_lengths_certify_enormous_sufficient_value_thresholds=False,
    later_generated_mean_and_repair_maps_automatically_identified=False,
    complete_curl_covariance_equals_principal_covariance=False,
    all_background_and_complete_residual_dependencies_closed=False,
    unconditional_pulse_family_or_blowup_proved=False,formal_certificate=False)
REMAINING=['consistently_join_canonical_primary_and_rebuilt_roots_to_current_generated_complete_covariance_fresh_mean_and_repair_maps',
    'complete_uniform_physical_scale_viscous_cutoff_interaction_residual_with_axial_velocity_and_full_pressure_memory',
    'retained_core_background_and_complete_construction_analytical_dependencies_before_unconditional_claims']


def protocol():
    return dict(parent_commit='817edc46cb8ba64db62f3e7c0287c45522ef5eab',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        canonical='equations 7.8 and Lemma 7.4: Br=U Mr(sref); zr(0)=P(0)*eplus, distinct from Ba entrance',
        fixed='same original phase, carriers, positive tilt, target, r0=1/6057072, width C^-20, cutoffs, gauge, unit viscosity; no historical pulse overwrite',
        canonical_column='CEr=CEref; D0r=Ds+Sref*Da+cmax*(1+Sref)*(Da+16CEr/(7mu)); D1=u*(1+cmax)',
        pulse='B0=2W0 max(1,cmax W0); Tstar=B0 exp(Cd+max(CEa,CEr)); CDelta=B0 Ds exp(Cd+CEa); |tr-ta|<=CDelta P/L',
        pressure='|pir-pia|<=Pstar CDelta P/(kL); full nprime and normal-source terms remain in forced inverse',
        Haar='s0=r0^2/sqrt(L); |DeltaH_sigma|<=s0 CH/L; CH=2Jmax Tstar CDelta sqrt(pi/c); angular half and actual cutoffs retained',
        roots='theta=12CH/(ja_minus*min(Amin,u)*L); kappa=24ja_plus/(gamma*ja_minus); theta<=1/(2kappa+1); relative squared-amplitude error R=kappa*theta/(1-theta)<=1/2',
        correction='yr=Hr^-1 T, ar=sqrt(yr); stale roots leave DeltaH*ya; equality Hr*yr=T is principal covariance only',
        domain='common normal, both frame Riccati and cone, and root-comparison VALUE thresholds; no jet-dependent band or giant numeric acceptance claim',
        memory='full moving-normal pressure; background Pi1/Pi2, axial velocity and whole pressure moment; Gaussian cutoff tails stay additive',
        scope='conditional canonical primary/root reconstruction and quantitative comparison, not automatic nonlinear map identification or complete residual')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'frame_column_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('canonical_current_frame_parent',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old


def history():return parent_module().parent_module().parent_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])


def sources():
    files=[HERE/n for n in ('README.md','audit.py','canonical_bounds.py','canonical_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def reviews():return deepcopy(dict(explicit_canonical_and_comparison_budgets=budget_review(),
    exact_principal_matrix_root_stability=matrix_review(),common_value_acceptance=acceptance_review(),
    independent_physical_frame_pressure_and_Haar_controls=physical_review()))


def compare_records(a,b,path='root'):
    if type(a) is not type(b):raise ValueError('Record type changed at '+path)
    if isinstance(a,dict):
        if a.keys()!=b.keys():raise ValueError('Record keys changed at '+path)
        for k in a:compare_records(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if len(a)!=len(b):raise ValueError('Record length changed at '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare_records(x,y,path+'.'+str(i))
    elif isinstance(a,float):
        tight=any(x in path for x in ('independent_frame_physical_error','normal_constraint','complete_pressure_balance_error','rebuilt_principal_covariance_error'))
        tol=2e-11 if tight else 2e-12
        if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=2e-6,abs_tol=tol):raise ValueError('Numerical diagnostic changed at '+path)
    elif a!=b:raise ValueError('Record changed at '+path)


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Canonical scope changed')
    compare_records(protocol(),r.get('protocol'),'protocol');compare_records(FLAGS,r.get('flags'),'flags')
    compare_records(REMAINING,r.get('remaining_obligations'),'limits')
    expected=reviews();compare_records(expected,r.get('reviews'),'reviews')
    g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in expected.items()};g['historical_recursive_source_bytes']=True
    compare_records(g,r.get('gates'),'gates')
    if not all(g.values()):raise ValueError('Canonical gate failed')
    compare_records(sources(),r.get('source_hashes'),'sources');compare_records(history(),r.get('inherited_provenance'),'history')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if a.verify_record:
        old=json.loads(a.verify_record.read_text());validate(old);compare_records(r,old)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Canonical principal family reconstructed conditionally; nonlinear map and full residual join remain.')
if __name__=='__main__':main()
