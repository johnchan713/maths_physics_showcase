#!/usr/bin/env python3
"""Current leading/frame floors, rounding, exact inverse and conditional scope."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from frame_floors import global_review,identity_review,logarithmic_review
from frame_checks import frame_review,carrier_review,pressure_and_domain_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-leading-frame-floors-from-reviewed-core-and-same-root-regional-profile'
DEPENDENCIES=('primary_covariance_foundation','global_pressure_moment_foundation','regional_cone_join','core_growth_activation','physical_pulse_foundation')
FLAGS=dict(uniform_leading_spectral_floors_from_current_reviewed_profile=True,
    specified_frame_left_inverse_bound_sqrt_two_G_squared=True,
    closed_outer_swirl_positive_despite_flat_zero_stress=True,
    angular_component_uniform_floor_not_assumed=True,
    axial_carrier_floor_keeps_both_signs_and_original_rounding=True,
    single_value_domain_for_every_fixed_frame_jet=True,
    current_primary_roots_and_sourced_inverse_keep_full_pressure=True,
    historical_evidence_unchanged_and_controls_not_selected_huge_profile=True,
    independently_formalized_infinite_core_contraction=False,
    actual_first_giant_band_numerically_evaluated=False,
    angular_frequency_implies_uniform_angular_component_floor=False,
    every_generated_complete_covariance_mean_and_residual_input_closed=False,
    unconditional_pulse_family_or_blowup_proved=False,formal_certificate=False)
REMAINING=['independent_review_of_retained_infinite_core_and_current_profile_analytical_maps',
    'complete_generated_covariance_and_fresh_radial_axial_mean_targets_with_Pi1_Pi2',
    'full_viscous_cutoff_interaction_residual_and_single_field_realization_before_unconditional_claims']


def protocol():
    return dict(parent_commit='889f2e5ba0837a3ea90abdacd7ae37d92f366b3d',
        fixed_choices='same Md64,h,N,matching lambda>0,unit viscosity,original tilt,carriers,primary roots,r0,frequency,width C^-20,supports,gauge',
        regional='alpha=1/(2G^3); nu=1/(2G^4); Ab=Bb=4Aenv; Vb=Ab+Bb^2/alpha',
        swirl='fstar=exp(-2Lambda*M)/(4C)*exp[-Ab*log(Xb/Xa)/2]; F<=2exp(2Lambda*M)/C',
        spectral='gstar^2=fstar^2*alpha*(2+nu); nstar^2=alpha/Vb; lambda_star^2=2alpha*fstar^2*nu/Vb; cstar^2=nu/2',
        chart='leading F,g,lambda multiply by s^(-1-h); N,K,c0 unchanged; 1/4<=s<=4',
        frame='specified Bleft=M^-1[er^T;Na^T]; norm2<=nu^-1/2=sqrt(2)*G^2',
        rounding='nearest NONZERO angular integer costs 1/k, not 1/(2k); no uniform |p/R| floor',
        axial='Ls>=2u/(gstar*nstar) gives pz<=-Bs*nstar/2; moving normal retains epsilon*v term',
        support='one common value domain before fixed jets; no derivative-dependent start; no new support or width',
        scope='leading/frame floors derived from reviewed current exact-core/regional maps, not all nonlinear inputs or blowup')


@lru_cache(None)
def parent_module():
    p=HERE.parent/DEPENDENCIES[0];before=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('leading_frame_current_parent',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=before


@lru_cache(None)
def parent_review():
    before=sys.path[:]
    try:
        m=parent_module();r=m.generate();m.validate(r)
        frozen=json.loads((HERE.parent/DEPENDENCIES[0]/'evidence.json').read_text());m.validate(frozen);m.compare_records(r,frozen)
    finally:sys.path[:]=before
    history=dict(r['inherited_provenance'])
    for name in DEPENDENCIES:
        record=json.loads((HERE.parent/name/'evidence.json').read_text())
        for path,digest in record['source_hashes'].items():history[name+':direct:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return r,history


def sources():
    paths=[HERE/n for n in ('README.md','audit.py','frame_floors.py','frame_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/n/'README.md' for n in ('core_growth_activation','regional_cone_join','heat_exterior_audit','outer_edge_stability','repair_state_bounds','loop_modulation_bounds')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(current_regional_and_logarithmic_swirl_certificate=global_review(),
    exact_leading_shear_and_squared_spectral_floors=identity_review(),logarithmic_and_open_domain_obstructions=logarithmic_review(),
    original_moving_frame_and_ambient_left_inverse=frame_review(),nonzero_rounding_and_axial_carrier=carrier_review(),
    full_pressure_and_common_domain_limits=pressure_and_domain_review()))


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    parent,history=parent_review();r=reviews();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g.update(current_primary_parent_reproduced=all(parent['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,inherited_provenance=history,
        source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Leading/frame scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Core or nonlinear obligations erased')
    expected=reviews();keys=set(expected)|{'current_primary_parent_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('Frame gate failed')
    compare_records(expected,r.get('reviews',{}),'reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if r.get('inherited_provenance')!=parent_review()[1] or not all(r['inherited_provenance'].values()):raise ValueError('History changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if a.verify_record:
        old=json.loads(a.verify_record.read_text());validate(old);compare_records(r,old)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Current leading/frame floors supplied from reviewed maps; generated complete covariance/mean/residual inputs remain open.')
if __name__=='__main__':main()
