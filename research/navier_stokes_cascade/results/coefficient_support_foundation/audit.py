#!/usr/bin/env python3
"""Current coefficient-support foundation with sharp linear jets and full pressure."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from support_bounds import HERE,load,map_review,jet_review,support_review
from conservative_checks import conservative_review,plateau_review
from linear_checks import numerical_review

PROJECT=HERE.parent.parent
STATUS='current-linear-coefficient-jets-reserved-velocity-and-exceptional-pressure-supports-joined-relative-to-leading-inputs'
FLAGS=dict(
    actual_positive_order_linear_map_and_continuum_inverse_800_supplied=True,
    all_fixed_eta_jet_normalization_costs_explicit=True,
    current_same_root_global_moment_inputs_replayed=True,
    every_finite_coefficient_support_and_reserved_mean_velocity_joined=True,
    all_order_conservative_totals_include_previous_order_axial_viscosity=True,
    first_two_pressure_tails_and_nonzero_leading_radial_memory_retained=True,
    exact_rational_checks_distinguished_from_manufactured_quadrature=True,
    historical_evidence_parameters_width_frequency_and_gauge_preserved=True,
    coefficient_norms_uniform_in_order_h_or_amplitude=False,
    reserved_mean_pressure_equal_to_leading_pressure=False,
    actual_selected_huge_profile_numerically_materialized=False,
    all_retained_leading_core_and_finite_pulse_lemmas_discharged=False,
    unconditional_actual_pulse_bounds_verified=False,
    unconditional_complete_covariance_mean_residual_or_blowup_verified=False,
    formal_certificate=False)
REMAINING=[
    'discharge_remaining_uniform_actual_pulse_existence_inverse_gap_size_derivative_and_support_lemmas_across_all_physical_scales',
    'retain_leading_radial_velocity_full_wave_axial_velocity_and_Pi1_Pi2_in_complete_covariance_mean_corrections',
    'prove_complete_residual_gain_including_viscosity_cutoffs_all_interactions_and_retained_core_dependencies_before_unconditional_claims']
DEPENDENCIES=('global_pressure_moment_foundation','coefficient_support_audit','weighted_stress_foundation')

def protocol():
    return dict(parent_commit='135c9da9a41f070ce0d583ba2a70d9a8fe21d490',
        fixed_choices='same Md64,lambda=exp(-4T)>0,h=lambda^2,N=1+floor(H^32),width C^-20,root branches,gauge and supports',
        ordinary_rows=['M','I','Cp','J','S'],
        positive_order_map='B*c_n=-z_n; no quadratic current-order term; ||B^-1||_infinity<800',
        amplitude='e(eta)=e_star/(1+eta^2), e_star>0; factorial reciprocal Ck costs 2,4,5 divided by e_star',
        jet_bound='||c_n||_k<=800 max(mu1/rho,(mu1/rho+K*Wk*mu4)/lambda,K*mu2,Wk*mu5/rho,Wk*mu3)',
        geometry='all I3 bumps end at y3=2.55; Xvel=rho3 exp(3); I4 has 6<y3<11',
        support_rows=['n1: velocity Xvel, pressure Xv, stress Xb','n2: velocity Xvel, pressure Xv, stress Xv','n>=3: velocity pressure stress Xvel'],
        mean_pressure=['Pi1=K1(eta)-mu^2/(4X)','Pi2=K2(eta)','Pi_n=0 for n>=3'],
        signed_totals='current moments plus shifted Z_(a-D)Z_a and Z_(g-D)Z_g of previous-order angular and axial moments',
        checks='48 independent exact-rational integrated identities; 8 non-small manufactured ordinary-integrand repairs and second jets',
        scope='coefficient/support foundation relative to retained finalized leading inputs; no uniform-in-order norms, unconditional pulse/residual/blowup or zero reserved pressure')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'global_pressure_moment_foundation';sys.path.insert(0,str(p))
    try:return load('coefficient_support_global_parent',p/'audit.py')
    finally:sys.path.pop(0)

@lru_cache(None)
def parent_review():
    module=parent_module();generated=module.generate();module.validate(generated)
    frozen=json.loads((HERE.parent/'global_pressure_moment_foundation/evidence.json').read_text())
    module.validate(frozen);module.compare_records(generated,frozen)
    history=dict(generated['inherited_provenance'])
    for p,d in frozen['source_hashes'].items():history['global_pressure_moment_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    for name in DEPENDENCIES[1:]:
        record=json.loads((HERE.parent/name/'evidence.json').read_text())
        for p,d in record['source_hashes'].items():history[name+':'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    return generated,history

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','support_bounds.py','conservative_checks.py','linear_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/p for p in ('annular_first_correction/extension.py','annular_first_correction/fixtures.py','background_summation/hierarchy.py','matching_root_foundation/matching_bounds.py')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():
    return deepcopy(dict(actual_linear_map=map_review(),all_fixed_normalized_eta_jets=jet_review(),
        common_and_exceptional_supports=support_review(),exact_conservative_all_order_identities=conservative_review(),
        shifted_radial_pressure_plateau=plateau_review(),actual_linear_integrand_diagnostics=numerical_review()))

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    parent,history=parent_review();r=reviews();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g.update(current_global_moment_parent_reproduced=all(parent['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
                inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Coefficient/support scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Remaining pulse/mean/residual obligation erased')
    if set(r.get('gates',{}))!=set(reviews())|{'current_global_moment_parent_reproduced','historical_bytes_preserved'} or not all(r['gates'].values()):raise ValueError('Support gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if r.get('inherited_provenance')!=parent_review()[1] or not all(r['inherited_provenance'].values()):raise ValueError('History changed')

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--verify-record',type=Path);args=parser.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if args.verify_record:
        f=json.loads(args.verify_record.read_text());validate(f);compare_records(r,f)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Finite coefficient supports and reserved velocity joined; Pi1/Pi2, uniform physical pulses and full residual limits retained.')
if __name__=='__main__':main()
