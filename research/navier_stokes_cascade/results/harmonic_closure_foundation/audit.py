#!/usr/bin/env python3
"""Explicit harmonic-budget audit. Provenance checks are separate from proofs."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from mode_bounds import convolution_review,inverse_review,linear_review,pressure_review,cycle_review
from mode_checks import complete_review,obstruction_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='complete-finite-harmonic-closure-with-explicit-mode-derivative-budget'
DEPENDENCIES=('leading_frame_foundation','sourced_inverse_foundation','primary_covariance_foundation',
              'compact_mean_foundation','regular_source_foundation','finite_induction_foundation')
FLAGS=dict(complete_finite_harmonic_convolution_and_covariance_count_independent=True,
    complete_transport_keeps_required_harmonic_ratio_derivative=True,
    sourced_velocity_two_index_and_pressure_one_index_smoothing_conditional=True,
    full_path_envelope_normalized_forcing_required_for_inverse=True,
    regular_linear_remainder_keeps_viscosity_and_gains_one_index=True,
    full_pressure_memory_current_axial_means_and_moving_moment_retained=True,
    same_common_band_original_geometry_and_exact_additive_tail_split=True,
    frozen_history_verified_by_hashes_separate_from_fresh_dependency_replays=True,
    exact_local_controls_not_actual_giant_profile_evaluation=True,
    generic_uniform_no_loss_harmonic_transport=False,
    actual_infinite_nonlinear_weighted_mode_sums_established=False,
    stage_uniform_constants_or_infinite_fourier_inverse_proved=False,
    every_retained_exact_core_and_construction_dependency_discharged=False,
    unconditional_complete_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['dependency_explicit_join_of_current_core_leading_primary_inverse_and_generated_complete_covariance_mean_maps',
    'complete_viscous_cutoff_interaction_residual_and_single_field_arguments_with_every_retained_dependency_discharged',
    'infinite_weighted_mode_sums_require_separate_uniform_majorants_if_an_infinite_Fourier_inverse_is_attempted']


def protocol():
    return dict(parent_commit='0c1399dca354f460b80c2b22ef80f63ea285fc95',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        fixed_choices='Md64; original h,N,lambda>0,unit viscosity,tilt,roots,carriers,width C^-20,supports,gauge',
        norm='N_sj=sum_m (1+abs(m))^s coefficient_jet_norm_j(a_m)',
        inverse='full-path envelope-normalized N^W: velocity gains 2 indices; potential 3 with k^-1; complete pressure 1 with k^-1; full entrance compatibility',
        transport='same-phase complete divergence removes k but leaves n/m; one angular derivative is necessary for generic uniform bounds',
        regular_linear='gains one harmonic index and epsilon^(1/2-kappa), with complete viscosity and pressure',
        covariance='all six symmetric entries and every matched negative-mode pair; old-new and new squares retained',
        memory='full Pi1,Pi2; current axial velocity and pressure; whole derivative Jz+c_rho*P; actual A0,A1 velocities',
        support='same accepted 4Q<qstar charts, full endpoint extension; exact Gaussian linear tails additive never forward-inverted',
        scope='finite count-independent estimates relative to reviewed maps; no infinite nonlinear summability or unconditional proof')


def verify_bindings(project,records):
    """Recursively verify frozen byte bindings, without running old diagnostics."""
    project=Path(project).resolve();todo=list(records);seen=set();out={}
    while todo:
        path=Path(todo.pop()).resolve()
        if path in seen:continue
        if not path.is_relative_to(project):raise ValueError('Record outside project')
        seen.add(path);record=json.loads(path.read_text())
        hashes=record.get('source_hashes',{})
        if not isinstance(hashes,dict):raise ValueError('Malformed historical bindings')
        for rel,digest in hashes.items():
            source=(project/rel).resolve()
            if not source.is_relative_to(project):raise ValueError('Source outside project')
            if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest()!=digest:raise ValueError('Historical source changed: '+rel)
            out[str(path.relative_to(project))+':'+rel]=True
            if source.name=='evidence.json':todo.append(source)
    return dict(records=len(seen),bindings=out)


def history():return verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])


def sources():
    paths=[HERE/n for n in ('README.md','audit.py','mode_bounds.py','mode_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(absolute_sum_convolution=convolution_review(),uniform_index_smoothing=inverse_review(),
    complete_regular_linear_budget=linear_review(),full_normal_pressure_limit=pressure_review(),finite_cycle_scope=cycle_review(),
    independent_complete_curl_transport_and_covariance=complete_review(),necessary_real_harmonic_loss=obstruction_review()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={key:bool(value['checks']) and all(value['checks'].values()) for key,value in r.items()}
    g['frozen_recursive_source_bindings']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):
    if type(a) is not type(b):raise ValueError('Record type changed at '+path)
    if isinstance(a,dict):
        if a.keys()!=b.keys():raise ValueError('Record keys changed at '+path)
        for k in a:compare_records(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if len(a)!=len(b):raise ValueError('Record length changed at '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare_records(x,y,path+'.'+str(i))
    elif a!=b:raise ValueError('Record changed at '+path)


def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Status changed')
    compare_records(protocol(),r.get('protocol'),'protocol');compare_records(FLAGS,r.get('flags'),'flags')
    compare_records(REMAINING,r.get('remaining_obligations'),'limits')
    expected=reviews();compare_records(expected,r.get('reviews'),'reviews')
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in expected.items()};gates['frozen_recursive_source_bindings']=True
    compare_records(gates,r.get('gates'),'gates')
    if not all(gates.values()):raise ValueError('Failed harmonic gate')
    compare_records(sources(),r.get('source_hashes'),'sources');compare_records(history(),r.get('inherited_provenance'),'history')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if a.verify_record:
        old=json.loads(a.verify_record.read_text());validate(old);compare_records(r,old)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Exact controls and frozen-byte checks; conditional analytical estimates, no unconditional proof.')
if __name__=='__main__':main()
