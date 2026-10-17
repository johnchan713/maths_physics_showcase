#!/usr/bin/env python3
"""Compact mean foundation with sharp fixed-lambda conditioning."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from mean_maps import reviews

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='compact-mean-inverse-foundation-sharp-axial-conditioning-relative-to-core-and-primary-inputs'
FLAGS=dict(all_fixed_order_compact_and_temporal_mean_inverse_estimates_derived=True,
    same_starting_domain_and_exact_physical_chart_intertwining=True,
    angular_five_map_nonsingular_at_lambda_zero=True,
    axial_one_over_lambda_loss_and_zero_limit_obstruction_derived=True,
    actual_positive_lambda_not_retuned=True,
    complete_covariance_axial_velocity_radial_flow_and_pressure_moments_retained=True,
    actual_background_Pi1_Pi2_and_flat_A0_A1_not_discarded=True,
    five_map_recomputed_nonlinear_defects_not_identically_zero=True,
    exact_core_and_remaining_background_pulse_lemmas_are_dependencies=True,
    historical_sources_evidence_supports_and_parameters_unchanged=True,
    exact_manufactured_checks_not_continuum_suprema=True,
    all_actual_generated_source_endpoint_and_coefficient_bounds_verified=False,
    full_later_complex_covariance_iteration_independently_closed=False,
    unconditional_all_pulse_family_bounds_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    actual_giant_parameters_numerically_evaluated=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['independently_verify_every_generated_nonlinear_harmonic_source_class_envelope_and_full_endpoint_extension',
    'independently_join_complete_later_complex_covariance_and_recomputed_mean_iterations',
    'derive_full_residual_with_viscosity_all_cutoffs_old_new_interactions_before_unconditional_claims']
DEPENDENCIES=('physical_pulse_foundation','compact_mean_corrections','complete_covariance_initial_means',
    'physical_wave_interface','coefficient_support_audit')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='372574a3e9f08bba9d31c351ce04298cd709a0ab',
        unchanged_parameters='Md64; unit viscosity; fixed lambda>0,h,N,t1=C^-20; original geometry and primary carriers',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        compact_remainder='norm_m(A_e f)<=C_mp*M^(-p)*norm_(m+2p+3)(f); zero radial Haar moment required',
        temporal_inverse='zero Haar input; norm_m(u)<=C_m*S*norm_(m+4)(f)',
        exact_floor_relation='M=c_i^rho_g*epsilon^(-kappa); M^(-1)<(Tg*S)^rho_g*epsilon^kappa',
        axial_scaled_identity='(1-q_lambda)*s1=Jtheta_q/(a*mu_(1-2lambda)*g); s0=-g*s1',
        physical_reserved_swirl='a(eta)*s^(lambda-h)*R^(-1-2lambda)',
        pressure_moment='int R*bar(p_m)=-int R^2*bar(g_r)/2+c_rho*P; differentiate WHOLE c_rho*P',
        scope='compact mean foundation and sharp conditioning, relative to core and primary inputs; all generated sources and complete residual remain open')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'physical_pulse_foundation';sys.path.insert(0,str(p))
    try:return load('compact_mean_foundation_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'physical_pulse_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    inherited=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():
        inherited['physical_pulse_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,inherited


def sources():
    paths=[HERE/n for n in ('README.md','mean_maps.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'compact_mean_corrections'/n for n in ('operators.py','mean_system.py','README.md')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=deepcopy(reviews())
    gates={k:bool(v) and all(v.values()) for k,v in r.items()}
    gates.update(actual_primary_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Mean foundation scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Source/full residual obligations erased')
    expected=deepcopy(reviews());keys=set(expected)|{'actual_primary_parent_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('A mean foundation gate failed')
    compare_records(expected,r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=874 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    args=p.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if args.verify_record:
        f=json.loads(args.verify_record.read_text());validate(f);compare_records(r,f)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Compact mean foundation only; all generated sources, later complex iteration, full residual and blow-up remain open.')


if __name__=='__main__':main()
