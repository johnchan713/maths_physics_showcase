#!/usr/bin/env python3
"""Join current same-root pressure, global moments and exterior stresses."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from identities import HERE,reviews as exact_reviews

PROJECT=HERE.parent.parent
STATUS='current-same-root-global-pressure-moments-and-leading-exterior-stresses-joined-relative-to-retained-inputs'
FLAGS=dict(
    transformed_matching_rows_proved_equivalent_to_original_five_moments=True,
    current_same_roots_and_positive_lambda_retained=True,
    complete_core_backward_pressure_same_gauge_trace_discharged=True,
    all_five_global_leading_moment_identities_joined=True,
    zero_defect_leading_exterior_signed_stress_identification_joined=True,
    inner_compact_outer_leading_value_cone_join_replayed_for_same_roots=True,
    full_pressure_memory_and_eta_derivative_retained=True,
    historical_frozen_evidence_preserved_and_replayed=True,
    retained_core_and_compact_input_estimates_independently_discharged=False,
    coefficient_and_reserved_mean_supports_verified=False,
    unconditional_uniform_actual_pulse_estimates_verified=False,
    complete_wave_covariance_mean_corrections_verified=False,
    unconditional_complete_residual_or_blowup_verified=False,
    formal_certificate=False)
REMAINING=[
    'discharge_remaining_finite_background_coefficient_and_reserved_mean_support_lemmas',
    'discharge_uniform_actual_pulse_existence_inverse_gap_size_derivative_and_support_estimates_at_every_physical_scale',
    'complete_wave_covariance_mean_corrections_including_axial_velocity_and_pressure',
    'prove_complete_residual_gain_with_viscosity_cutoffs_and_all_interactions_before_unconditional_claims']
DEPENDENCIES=('matching_root_foundation','global_moment_stress','regional_cone_join')


def protocol():
    return dict(parent_commit='29b1f5f281650f5fae082d5d4a90a1f0d66e9ae3',
        parameters='same Md64,T=exp(64)+10,lambda=exp(-4T),h=lambda^2,Tw=240T,Tf=1000,co=.001,width C^-20',
        rows='original global defects (M,R,J,S,Cp); R is the renormalized heat moment with the unchanged final power tail',
        annular_transform='(M,J-4eta I,I,-(S-8eta M),Cp), with nonzero positive normalizers omitted only for zero-row equivalence',
        compact_transform='(M,(M-J)/lambda,I,-S,Cp), with retained 0<lambda<=.01',
        pressure='Pi_axis=-Cp_id(infinity); same Pi_axis+Cp gauge at every join; no reset; Pi=-tail Cp',
        exterior='zero five functional defects imply zero V0 and signed weighted stress defects; full pressure X term and eta derivative retained',
        cone='same accepted annular and compact roots, unchanged N and C^-20, so the frozen regional leading-value cone join applies',
        scope='current same-root global leading pressure/moment/stress join relative to retained core and compact-input estimates; finite-background/pulse and complete covariance/residual obligations remain')


def _load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    sys.path.insert(0,str(path.parent))
    try:spec.loader.exec_module(module)
    finally:sys.path.pop(0)
    return module


@lru_cache(None)
def dependency_review(name):
    folder=HERE.parent/name;module=_load('global_pressure_dependency_'+name,folder/'audit.py')
    generated=module.generate();module.validate(generated)
    frozen=json.loads((folder/'evidence.json').read_text());module.validate(frozen);module.compare_records(generated,frozen)
    history={name+':inherited:'+k:v for k,v in frozen['inherited_provenance'].items()}
    for p,d in frozen['source_hashes'].items():
        history[name+':source:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    return generated,history


def source_hashes():
    paths=[HERE/n for n in ('README.md','identities.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'global_moment_stress/ledger.py',HERE.parent/'regional_cone_join/cones.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generated_reviews():
    r=exact_reviews()
    matching,_=dependency_review('matching_root_foundation')
    global_mom,_=dependency_review('global_moment_stress')
    regional,_=dependency_review('regional_cone_join')
    r['current_chain_scope']=dict(checks=dict(
        current_matching_same_roots_reproduced=all(matching['gates'].values()),
        historical_global_moment_derivation_reproduced=all(global_mom['gates'].values()),
        regional_leading_value_cone_join_reproduced=all(regional['gates'].values()),
        current_chain_has_all_five_pressure_memory_rows=matching['flags']['all_five_moments_and_complete_pressure_memory_retained'],
        same_finite_frequency_cone_joined=regional['flags']['same_finite_frequency_leading_value_cone_joined'],
        compact_input_cone_hypotheses_discharged=regional['flags']['compact_input_cone_hypotheses_discharged'],
        downstream_pulse_and_residual_claims_remain_conditional=True))
    return r


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=generated_reviews();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    history={}
    for name in DEPENDENCIES:history.update(dependency_review(name)[1])
    g['historical_bytes_preserved']=bool(history) and all(history.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
                inherited_provenance=history,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):
    if type(a) is not type(b):raise ValueError('Changed type at '+path)
    if isinstance(a,dict):
        if set(a)!=set(b):raise ValueError('Changed keys at '+path)
        for k in a:compare_records(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if len(a)!=len(b):raise ValueError('Changed list at '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare_records(x,y,path+'.'+str(i))
    elif isinstance(a,float):
        if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-65):raise ValueError('Changed float at '+path)
    elif a!=b:raise ValueError('Changed value at '+path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Remaining obligation erased')
    if set(r.get('reviews',{}))!=set(generated_reviews()):raise ValueError('Review removed')
    if set(r.get('gates',{}))!=set(generated_reviews())|{'historical_bytes_preserved'} or not all(r['gates'].values()):raise ValueError('Gate failed')
    compare_records(generated_reviews(),r['reviews'],'root.reviews')
    if r.get('source_hashes')!=source_hashes():raise ValueError('Source bytes changed')
    expected={}
    for name in DEPENDENCIES:expected.update(dependency_review(name)[1])
    if r.get('inherited_provenance')!=expected or not all(expected.values()):raise ValueError('Historical provenance changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Current same-root leading pressure/moment/stress join supplied; finite-pulse and complete-residual obligations remain.')
if __name__=='__main__':main()
