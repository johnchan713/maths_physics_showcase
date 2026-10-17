#!/usr/bin/env python3
"""First actual sourced pulse, complete covariance and recomputed pressure."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from sourced_bounds import review
from sourced_checks import reviews as exact_reviews

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='first-sourced-physical-pulse-complex-covariance-and-current-pressure-foundation-relative-to-core-and-finite-lemmas'
FLAGS=dict(unchanged_signed_mean_prefix_and_same_first_inverse_source_preserved=True,
    actual_first_sourced_pulse_all_fixed_weighted_and_physical_jets_joined=True,
    all_accepted_scales_same_band_domain_closed_supports_and_endpoint_extension=True,
    complex_source_normal_and_moving_normal_pressure_terms_retained=True,
    all_six_complex_covariances_and_second_harmonic_square_retained=True,
    full_pressure_and_three_current_defects_recomputed_with_nonzero_old_means=True,
    full_mean_order_D_Haar_mean_and_defect_order_C_distinguished=True,
    original_background_Pi1_Pi2_and_whole_moving_bump_product_retained=True,
    one_finite_regular_nonzero_gain_one_half_minus_kappa_derived=True,
    actual_flat_velocities_in_all_products_and_Gaussian_tails_not_inverse_sources=True,
    exact_core_and_finite_background_pulse_lemmas_remain_dependencies=True,
    historical_evidence_parameters_geometry_width_and_original_roots_unchanged=True,
    manufactured_diagnostics_not_actual_giant_field_or_continuum_suprema=True,
    initial_real_covariance_shortcut_applied_to_sourced_complex_pulse=False,
    generic_covariance_increment_claimed_positive_semidefinite=False,
    old_wave_error_promoted_to_new_residual_order=False,
    all_later_complex_targets_and_full_induction_verified=False,
    unconditional_all_pulse_family_bounds_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    actual_giant_profile_or_first_band_numerically_evaluated=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['join_current_generic_complex_signed_target_order_C_with_old_error_rho_B',
    'relax_current_fast_and_slow_means_recomputing_full_pressure_and_defects_after_each',
    'independently_close_later_source_families_and_complete_residual_before_unconditional_claims']
DEPENDENCIES=('regular_source_foundation','physical_pulse_foundation','compact_mean_foundation',
    'real_signed_covariance_update','post_signed_mean_relaxation','first_sourced_wave_update',
    'complete_covariance_initial_means','complex_signed_covariance_update')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='9dbaa3cf768f3009f4f17d9b5c5b47e4df613a6b',
        unchanged_parameters='Md64; unit viscosity; fixed h,N,positive lambda,t1=C^-20; original radius roots supports carriers and operation order',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_order='B=1-kappa; same source after unchanged real signed and post-signed mean prefix',
        covariance_pressure_defect_order='C=B+1/2=3/2-kappa',
        regular_nonzero_and_full_mean_order='D=C-kappa=3/2-2*kappa',
        complete_old_wave_error='rho=B, not D; original W0_as is fixed',
        pressure='i*(n.Kt-nprime.t+n.f)/(k*m*|n|^2); pressure and potential order C',
        physical_powers='u:-A+hB; p:-2A+hC; potential:1/2-A+hC; subtract a*sx+b*st',
        exact_tail_split='psi*uncut_regular_remainder + (1-psi)*f + T_psi; tails additive not forward sources',
        scope='ONE finite sourced/covariance/pressure join relative to exact core and finite lemmas; later complex induction and unconditional full residual remain open')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'regular_source_foundation';sys.path.insert(0,str(p))
    try:return load('first_sourced_foundation_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'regular_source_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():history['regular_source_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,history


def sources():
    paths=[HERE/n for n in ('README.md','sourced_bounds.py','sourced_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'complete_covariance_initial_means/means.py',HERE.parent/'actual_band_pulse_bounds/bounds.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(same_prefix_first_sourced_and_physical_orders=review(),**exact_reviews()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews();gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(actual_source_parent_and_retained_chain_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('First sourced foundation scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Later complex/full residual obligations erased')
    keys=set(reviews())|{'actual_source_parent_and_retained_chain_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('A first sourced foundation gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=907 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);args=p.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if args.verify_record:
        f=json.loads(args.verify_record.read_text());validate(f);compare_records(r,f)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('One finite first sourced join; all later complex states, full residual and blow-up remain open.')


if __name__=='__main__':main()
