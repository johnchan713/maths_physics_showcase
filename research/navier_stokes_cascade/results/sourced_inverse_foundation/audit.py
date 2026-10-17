#!/usr/bin/env python3
"""Explicit sourced inverse, with source compatibility and pressure limits."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from inverse_bounds import analytical_review
from inverse_checks import principal_checks,numerical_checks,endpoint_checks

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='explicit-same-domain-sourced-inverse-and-high-harmonic-parameter-jets-relative-to-reviewed-frame-inputs'
DEPENDENCIES=('coefficient_support_foundation','physical_pulse_foundation')
FLAGS=dict(explicit_finite_sourced_ODE_existence_uniqueness_from_reviewed_frame=True,
    high_harmonic_velocity_parameter_jets_inverse_square_derived=True,
    entrance_compatible_pulse_jets_have_same_high_mode_gain=True,
    source_normal_pressure_only_inverse_first_harmonic_power_retained=True,
    same_accepted_domain_for_all_nonzero_integer_modes=True,
    actual_background_radial_flow_and_Pi1_Pi2_not_erased=True,
    generic_sources_do_not_get_automatic_pulse_jet_gain=True,
    uncut_exit_not_declared_compactly_supported=True,
    fixed_finite_harmonic_physical_transfer_not_infinite_sum=True,
    historical_evidence_and_parameters_unchanged=True,
    manufactured_controls_not_actual_huge_profile=True,
    all_leading_core_frame_primary_root_dependencies_discharged=False,
    all_actual_nonlinear_source_and_mean_inputs_independently_closed=False,
    unconditional_all_actual_pulse_family_verified=False,
    unconditional_full_residual_or_blowup_verified=False,
    infinite_harmonic_physical_summability_verified=False,
    formal_certificate=False)
REMAINING=[
    'independently_join_actual_generated_source_and_complete_covariance_mean_inputs_to_current_roots_pressure_and_supports',
    'retain_exact_core_actual_background_frame_and_original_positive_primary_root_dependencies',
    'verify_full_viscous_cutoff_interaction_residual_and_single_field_realization_before_unconditional_claims']


def protocol():
    return dict(parent_commit='93103d22e03b2996ad8cf18ef31f77d090c22896',
        fixed_choices='same Md64,h,N,positive lambda,original roots,r0,width C^-20,frequency,supports and accepted band domain',
        frame_inputs='lambda>0; d>=dstar>0; ||E||+|d-dref|<=H/L; all fixed coefficient jets polynomial in S',
        propagator='||Vm(v,w)||<=exp(H)*P(v)/P(w)*exp[-(m^2-1)dstar*(v-w)]',
        residence='ell_m=L for |m|=1; min(L,1/((m^2-1)dstar)) for |m|>=2',
        high_jet_recursion='ZI=C0[FI/dstar+sum binomial(I,J)(RJ/(3dstar)+dJ/dstar)Z_(I-J)]',
        pulse_jet_condition='all required entrance source jets vanish; full smooth temporal zero extension suffices',
        pressure='i*(n.Kt-nprime.t+n.f)/(k*m*|n|^2); direct normal source only 1/(k*m)',
        supports='full lifted path parameter supports; uncut exit may be nonzero; cut potential AND pressure before curl',
        physical_powers='-A+h*alpha-a*sx-b*st; -2A+h*(alpha+1/2)-a*sx-b*st; 1/2-A+h*(alpha+1/2)-a*sx-b*st',
        scope='explicit linear sourced inverse from retained frame inputs, not all nonlinear sources, infinite harmonics, unconditional pulses, residual or blowup')


@lru_cache(None)
def parent_module(name):
    if name not in DEPENDENCIES:raise ValueError('Unknown dependency')
    p=HERE.parent/name;before=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('sourced_inverse_parent_'+name,p/'audit.py')
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
        spec.loader.exec_module(module);return module
    finally:sys.path[:]=before


@lru_cache(None)
def parent_review():
    records={};history={}
    for name in DEPENDENCIES:
        before=sys.path[:]
        try:
            module=parent_module(name);record=module.generate();module.validate(record)
            frozen=json.loads((HERE.parent/name/'evidence.json').read_text())
            module.validate(frozen);module.compare_records(record,frozen)
        finally:sys.path[:]=before
        records[name]=record
        for key,value in record['inherited_provenance'].items():history[name+':'+key]=value
        for path,digest in frozen['source_hashes'].items():history[name+':direct:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return records,history


def sources():
    paths=[HERE/name for name in ('README.md','audit.py','inverse_bounds.py','inverse_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/name/'evidence.json' for name in DEPENDENCIES]
    paths += [HERE.parent/name for name in ('regular_source_foundation/source_bounds.py','physical_pulse_foundation/pulse_bounds.py','first_sourced_foundation/sourced_bounds.py')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():
    return deepcopy(dict(analytical_inverse_and_high_mode_jets=analytical_review(),
        exact_moving_principal_pressure_controls=principal_checks(),
        manufactured_noncommuting_ODE_and_scalar_parameter_controls=numerical_checks(),
        entrance_compatibility_and_temporal_support=endpoint_checks()))


def compare_records(a,b,path='root'):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b):raise ValueError('Record keys changed at '+path)
        for k,v in a.items():compare_records(v,b[k],path+'.'+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise ValueError('Record list changed at '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare_records(x,y,path+'.'+str(i))
    elif isinstance(a,float):
        if not isinstance(b,(int,float)) or isinstance(b,bool) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-10):raise ValueError('Numerical diagnostic changed at '+path)
    elif type(a) is not type(b) or a!=b:raise ValueError('Exact record changed at '+path)


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    parents,history=parent_review();r=reviews()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(current_coefficient_and_physical_parents_reproduced=all(all(p['gates'].values()) for p in parents.values()),
        historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Sourced inverse scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Remaining obligations erased')
    expected=reviews();keys=set(expected)|{'current_coefficient_and_physical_parents_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('Inverse gate failed')
    compare_records(expected,r.get('reviews',{}),'reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if r.get('inherited_provenance')!=parent_review()[1] or not all(r['inherited_provenance'].values()):raise ValueError('History changed')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--verify-record',type=Path);args=parser.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text());validate(frozen);compare_records(r,frozen)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Explicit sourced inverse from reviewed frame; full nonlinear inputs and unconditional residual/blowup remain open.')
if __name__=='__main__':main()
