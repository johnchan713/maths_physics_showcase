#!/usr/bin/env python3
"""Guarded exact value-domain certificate and unchanged-first-band limit."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from acceptance import load,dyadic_review,domain_review,logarithmic_review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='conditional-dyadic-value-domain-certificate-not-first-band-verification'
DEPENDENCIES=('core_pressure_identification','canonical_primary_foundation','frame_column_foundation',
    'leading_frame_foundation','geometric_pulse_replay','canonical_endpoint_join')
FLAGS=dict(exact_all_later_band_value_certificate_derived=True,
    canonical_direct_roots_separated_from_historical_relative_comparison=True,
    positive_spectral_floors_traced_to_current_same_profile=True,
    actual_finite_background_norm_dependence_retained=True,
    original_radius_Haar_cost_and_whole_pulse_rounding_retained=True,
    one_value_domain_for_every_fixed_jet=True,
    old_qstar_fields_supports_cutoffs_and_schedules_unchanged=True,
    complete_normal_source_and_Pi1_Pi2_pressure_memory_retained=True,
    parameter_fixtures_distinguished_from_actual_enclosures=True,
    actual_outward_norm_enclosures_numerically_verified=False,
    original_first_band_or_old_localized_force_value_acceptance_verified=False,
    failed_sufficient_comparison_gate_implies_nonexistent_canonical_roots=False,
    observation_tail_identified_with_old_complete_domain=False,
    new_localized_force_silently_substituted_for_old=False,
    all_retained_analytical_construction_premises_discharged=False,
    unconditional_pulse_complete_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['verify_actual_outward_value_norm_enclosures_and_existing_minimum_band_inequality_or_explicitly_propagate_a_distinct_certified_tail_theorem_through_localization',
    'discharge_actual_core_background_existence_inverse_edge_and_finite_contraction_premises_before_unconditional_claims']

def protocol():
    return dict(parent_commit='edf21536d1268a5d8061d663427128c361ec4970',
        certificate='Ls>=Lreq<=2^M; ell>=2^B; B=max(4,b(ell_old),b(512/h^2),ceil((M+b(1/(2r0)))/2))',
        exact='b(x)=ceil(log2(x)) by exact integer bit lengths; exp(x)<=2^(2x); sqrt(pi/c)<=2(1+1/c)',
        smallness='3 ell^4 exp(-h ell/4)<=3*8!*4^8/(h^8 ell^4)<1 for ell>=512/h^2',
        scale='Ls=2r0/ci>=2r0 ell^2, not log(1/Q); original radius and width C_setup^-20 retained',
        distinction='canonical own cone gives positive roots; relative comparison with historical roots is stronger and not necessary for existence',
        domain='qtail=min(qstar/64,2^(-E-2)), E=2^B, is observation only; old field and localized force never replaced',
        missing='actual outward norm enclosures and old minimum contributing band have not been supplied numerically; eventual acceptance does not verify old first band',
        pressure='complete moving-normal/source-normal pressure, axial means and Pi1/Pi2 memory unchanged',
        scope='conditional analytical value certificate; exact rational and mpmath fixtures not actual construction or continuum proof')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'core_pressure_identification';old=sys.path[:];sys.path.insert(0,str(p))
    try:return load('value_domain_pressure_parent',p/'audit.py')
    finally:sys.path[:]=old

@lru_cache(None)
def history():return parent_module().parent_module().provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','acceptance.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    paths += [HERE.parent/'canonical_primary_foundation/canonical_bounds.py',HERE.parent/'frame_column_foundation/frame_budget.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

@lru_cache(None)
def reviews():return dict(exact_dyadic_budgets=dyadic_review(),domain_and_gate_distinctions=domain_review(),independent_logarithmic_diagnostics=logarithmic_review())

def compare_records(a,b,path='root'):
    # All serialized new controls, including diagnostic booleans, reproduce exactly.
    if type(a)!=type(b):raise ValueError('Type changed at '+path)
    if isinstance(a,dict):
        if a.keys()!=b.keys():raise ValueError('Keys changed at '+path)
        for k in a:compare_records(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if len(a)!=len(b):raise ValueError('Length changed at '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare_records(x,y,path+str(i))
    elif a!=b:raise ValueError('Value changed at '+path)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    expected=generate();compare_records(expected,r)
    if not all(expected['gates'].values()):raise ValueError('Value-domain gate failed')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('All-later value certificate only; actual old-band enclosures, remaining analytical premises and unconditional proof open.')
if __name__=='__main__':main()
