#!/usr/bin/env python3
"""Explicit phase/frame/column constants with exact history and scope guards."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from frame_budget import budget_review,all_band_review
from frame_identity import identity_review,eigen_review,column_review,gauge_review
from phase_checks import phase_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='explicit-current-phase-frame-and-primary-column-comparison-budgets'
DEPENDENCIES=('gaussian_mass_foundation','leading_frame_foundation','primary_covariance_foundation',
              'physical_pulse_foundation','geometric_pulse_replay')
FLAGS=dict(original_full_phase_gives_explicit_primitive_normal_and_jet_errors=True,
    exact_two_coordinate_frame_generator_keeps_nprime_and_Bprime=True,
    specified_ambient_left_inverse_used_without_pseudoinverse_substitution=True,
    explicit_Euclidean_operator_CE_and_full_damping_Cd=True,
    complete_principal_column_ratio_keeps_secondary_frame_coordinate=True,
    explicit_D0_D1_and_logarithmic_Cdir_in_existing_cone_budget=True,
    same_accepted_value_domain_for_every_fixed_jet_and_later_source=True,
    original_radius_carriers_roots_width_supports_and_pressure_memory_retained=True,
    rational_and_high_precision_controls_not_continuum_certificates=True,
    manuscript_reference_mixing_distinguished_from_inherited_actual_mixing=True,
    physical_pulse_invariant_under_coordinate_change_with_transported_entrance=True,
    canonical_growing_entrances_in_both_frames_identified_as_same_pulse=False,
    huge_actual_background_seminorms_or_first_band_numerically_evaluated=False,
    arbitrary_existing_numeric_first_band_declared_accepted=False,
    positive_full_oscillating_velocity_lower_bound_claimed=False,
    every_generated_covariance_fresh_mean_and_full_residual_input_closed=False,
    retained_core_and_background_analytical_dependencies_discharged=False,
    unconditional_construction_or_blowup_proved=False,formal_certificate=False)
REMAINING=['resolve_manuscript_canonical_primary_entrance_versus_inherited_actual_mixing_before_identifying_their_exact_roots',
    'join_explicit_current_core_background_frame_primary_inverse_and_generated_complete_covariance_fresh_mean_maps',
    'complete_viscous_cutoff_interaction_residual_on_same_current_domain_and_realization_with_all_dependencies_accounted',
    'retained_exact_core_and_background_analytical_dependencies_before_unconditional_claims']


def protocol():
    return dict(parent_commit='255ae2e44292419235dd9b7e1a340f20e690dc39',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        fixed_choices='Md64; original h,N,positive tilt,lambda>0,roots,carriers,r0,width C_setup^-20,supports,gauge,unit viscosity',
        primitive='same phase; full 1/k rounding and epsilon*v axial term; Db,Dbox,Mz from actual reviewed background jets',
        exact_frame='H11=s*(a.g-sprime)/(1+s^2); H12=(2F*btheta-s*omega)/(1+s^2); H21=-(2F*btheta+b.g)+s*omega; H22=0',
        connection='Bleft*(An*B-Bprime)=M^-1*H*M-M^-1*Mprime; M^-1*Mprime=ssprime/[2(1+s^2)]*[[1,-1],[-1,1]]',
        frame_convention='sections 2-4 use inherited actual-sa mixing; manuscript (7.8) uses sref mixing; T=Ma^-1*Mr; zr(0)=T(0)^-1*za(0), not reset eplus',
        manuscript_budget='CEref=kMref*(DH+2*Fmax*Ds)+Sref*u; full time-dependent generator and forcing transform',
        operator='CE=kM*DH+S0*S1+lambda_plus*Ds; Cd=4*Dn*(2*Nref+Dn); Euclidean operator not entrywise norm',
        column='ttan/x=-s*a+c0*w*(1-q)/(1+q)*b; |q|<=1/8; D0 and D1 explicit; Cdir=D0+D1*Gaussian_moment_factor',
        acceptance='original normal floors, L>=8CE/mu, L>=25 and primary-cone value inequalities; no derivative-dependent band or numeric first-band claim',
        memory='full nprime and source-normal pressure; Pi1/Pi2; axial velocity and whole moving pressure moment; additive Gaussian tails never forward-inverted',
        scope='explicit frame/column comparison relative to current core/background jet bounds; complete covariance/mean/residual join and unconditional proof remain')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'gaussian_mass_foundation';before=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('frame_column_current_gaussian',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=before


def history():return parent_module().parent_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])


def sources():
    files=[HERE/n for n in ('README.md','audit.py','frame_budget.py','frame_identity.py','phase_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def reviews():return deepcopy(dict(explicit_current_phase_frame_column_budgets=budget_review(),same_value_band_domain=all_band_review(),
    independent_exact_projected_frame_identity=identity_review(),original_reference_diagonalization=eigen_review(),
    complete_principal_column_ratio=column_review(),manuscript_repository_frame_distinction=gauge_review(),
    original_full_phase_high_precision_controls=phase_review()))


def compare_records(a,b,path='root'):
    if type(a) is not type(b):raise ValueError('Record type changed at '+path)
    if isinstance(a,dict):
        if a.keys()!=b.keys():raise ValueError('Record keys changed at '+path)
        for k in a:compare_records(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if len(a)!=len(b):raise ValueError('Record length changed at '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare_records(x,y,path+'.'+str(i))
    elif isinstance(a,float):
        tolerance=1e-65 if path.endswith('maximum_direct_reduced_identity_error') else 2e-10
        if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=2e-8,abs_tol=tolerance):raise ValueError('Numerical control changed at '+path)
    elif a!=b:raise ValueError('Record changed at '+path)


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Frame/column scope changed')
    compare_records(protocol(),r.get('protocol'),'protocol');compare_records(FLAGS,r.get('flags'),'flags')
    compare_records(REMAINING,r.get('remaining_obligations'),'limits')
    expected=reviews();compare_records(expected,r.get('reviews'),'reviews')
    g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in expected.items()};g['historical_recursive_source_bytes']=True
    compare_records(g,r.get('gates'),'gates')
    if not all(g.values()):raise ValueError('Frame/column gate failed')
    compare_records(sources(),r.get('source_hashes'),'sources');compare_records(history(),r.get('inherited_provenance'),'history')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if a.verify_record:
        old=json.loads(a.verify_record.read_text());validate(old);compare_records(r,old)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Explicit conditional frame/column estimates; complete residual and unconditional proof remain open.')
if __name__=='__main__':main()
