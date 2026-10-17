#!/usr/bin/env python3
"""Local all-order leading-core foundation; global profile/pulse claims remain conditional."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from foundation import analytic_review,zero_free_review,profile_review,physical_review
from coordinate_checks import diagnostic_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='local-core-all-order-joint-analytic-reciprocal-and-Cartesian-scale-foundation'
FLAGS=dict(infinite_core_coefficient_solution_not_finite_iteration_is_input=True,
    explicit_all_order_joint_majorant_derived=True,
    common_zero_free_reciprocal_neighborhood_and_all_fixed_jets_derived=True,
    exact_core_pressure_datum_increment_and_parameter_derivatives_retained=True,
    Cartesian_axis_extension_and_anisotropic_all_scale_core_jet_bounds_derived=True,
    physical_exponent_updated_after_each_coordinate_derivative=True,
    actual_unchanged_leading_core_scope_is_Y_at_most_four=True,
    manufactured_diagnostics_are_not_the_enormous_core=True,
    reference_pressure_family_and_local_contraction_hypotheses_explicit=True,
    full_global_foundation_independently_verified=False,
    common_analytic_collar_after_flat_activation_independently_verified=False,
    all_modulated_repaired_leading_edge_bounds_independently_verified=False,
    unconditional_actual_pulse_bounds_verified=False,
    complete_viscous_residual_flatness_follows_from_core_stress_cancellation=False,
    constants_uniform_in_derivative_order_or_as_C_tends_to_infinity=False,
    singular_event_velocity_extension_verified=False,
    unconditional_blowup_verified=False,formal_certificate=False)
REMAINING=['extend_the_same_all_order_leading_certificate_through_actual_flat_activation_to_the_background_collar',
    'independently_audit_modulation_repair_global_pressure_moments_and_weighted_edges_from_defining_maps',
    'then_validate_all_order_background_pulse_inverse_and_physical_estimates_before_unconditional_claims']
DEPENDENCIES=('current_solution_comparison','core_contraction_audit','core_growth_activation','background_jet_foundation','paper_profile_audit')


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='19e184739d019e14c3e151a9f3a49dbe7a2d9f65',
        source='repository core_contraction_audit defining coefficient norm and core map; no final manuscript theorem is input',
        retained_pressure='exact positive measure on [0,1], mass<5, P>=16; full Pi0 plus p/Lambda',
        h_range='0<h<=1/100',rho='sigma/2048',same_Lambda_C_h_N_and_width=True,
        joint_majorant='norm*(r+k)!*20^(-r)*rho^(-k)*(1-|Y|/20-|zeta|/rho)^(-r-k-1)',
        majorant_domain='|Y|<=5, complex eta displacement<=rho/8',
        reciprocal_polydisc=['1/(100*S)','rho/(100*S)'],reciprocal_upper='5',
        seed_real_interval='0<=Y<=4.1',actual_unchanged_leading_interval='0<=Y<=4',
        physical_scale='0<q<=1; q=0 singular event excluded',physical_weight='(alpha_x+alpha_y)/2+(1/2-h)*alpha_z+j',
        diagnostic_h='1/200',diagnostic_Lambda=3,diagnostic_order=4,diagnostic_arithmetic='exact rational truncated multivariate Taylor algebra',
        scope='local core foundation only; actual flat activation collar, global leading profile, full background/pulses and blow-up remain unverified')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'current_solution_comparison';sys.path.insert(0,str(p))
    try:return load('core_foundation_current_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    frozen=json.loads((HERE.parent/'current_solution_comparison/evidence.json').read_text())
    m.validate(frozen);m.compare_records(r,frozen)
    inherited=dict(r['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['current_solution_comparison:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return r,inherited


@lru_cache(None)
def core_review():
    p=HERE.parent/'core_contraction_audit'
    saved={n:sys.modules.pop(n,None) for n in ('core','checks')};sys.path.insert(0,str(p))
    try:
        m=load('core_foundation_infinite_core',p/'audit.py');r=m.generate();m.validate(r)
        frozen=json.loads((p/'evidence.json').read_text());m.validate(frozen);m.compare_records(r,frozen)
    finally:
        sys.path.pop(0)
        for n,v in saved.items():
            sys.modules.pop(n,None)
            if v is not None:sys.modules[n]=v
    return dict(controls={k:{n:bool(v) for n,v in a['checks'].items()} for k,a in r['reviews'].items()},
        checks=dict(infinite_space_contraction_and_original_pressure_equations_reproduce=all(r['gates'].values()),
            core_source_bytes_unchanged=r['source_hashes']==frozen['source_hashes'],
            pressure_family_is_still_an_explicit_local_theorem_hypothesis=True,
            local_stress_cancellation_not_replaced_by_zero_complete_residual=True))


def sources():
    names=('README.md','foundation.py','coordinate_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')
    paths=[HERE/n for n in names]+[HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():
    return dict(all_index_joint_analytic_majorant=analytic_review(),uniform_zero_free_reciprocal=zero_free_review(),
        complete_core_profile_and_axis_quotient=profile_review(),all_scale_Cartesian_core_bounds=physical_review(),
        exact_independent_coordinate_and_failure_controls=diagnostic_review(),infinite_core_premises_reproduced=deepcopy(core_review()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    parent,inherited=parent_review();r=reviews()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates['current_conditional_chain_reproduced_without_promotion']=all(parent['gates'].values())
    gates['historical_source_bytes_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():raise ValueError('Local core scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Global foundation gate erased')
    expected=reviews();keys=set(expected)|{'current_conditional_chain_reproduced_without_promotion','historical_source_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('A core foundation gate failed')
    if r.get('source_hashes')!=sources():raise ValueError('Current source bytes changed')
    if len(r.get('inherited_provenance',{}))!=758 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(expected,r.get('reviews',{}),'root.reviews')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    args=p.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; historical frozen evidence is never overwritten')
    r=generate();validate(r)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text());validate(frozen);compare_records(r,frozen)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates and {len(r['inherited_provenance'])} inherited hashes pass")
    print('Local pressure-family core only; flat activation collar, global foundation, unconditional pulses and blow-up remain open.')


if __name__=='__main__':main()
