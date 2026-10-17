#!/usr/bin/env python3
"""Current-input complex signed correction followed by refreshed mean maps."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from complex_join_bounds import review
from complex_join_checks import reviews as exact_reviews

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-input-complex-signed-and-fresh-mean-foundation-relative-to-core-and-finite-lemmas'
FLAGS=dict(same_current_complex_signed_fast_slow_operations_and_original_roots=True,
    current_B_C_D_inputs_joined_to_all_fixed_signed_and_mean_estimates=True,
    all_three_generic_covariance_error_groups_retained_and_tie=True,
    full_updated_tensor_and_pressure_recomputed_after_each_operation=True,
    fresh_postfast_defects_used_by_same_fixed_positive_lambda_five_map=True,
    old_Haar_means_axial_velocity_pressure_and_actual_A0_A1_retained=True,
    whole_moving_pressure_product_and_both_Pi1_Pi2_tails_retained=True,
    fast_mean_wave_source_only_ties_existing_D_not_strictly_higher=True,
    full_tangential_mean_nu_and_fresh_defects_C_plus_one_joined=True,
    all_accepted_physical_scales_closed_supports_and_flat_tail_ledger_retained=True,
    exact_core_and_finite_defining_inverse_lemmas_remain_dependencies=True,
    historical_frozen_bytes_parameters_width_geometry_and_operation_order_unchanged=True,
    manufactured_checks_not_actual_giant_field_or_continuum_suprema=True,
    initial_real_phase_shortcut_applied_to_complex_state=False,
    whole_covariance_change_promoted_to_nu=False,
    stale_prefast_targets_suffice_for_final_defect_bound=False,
    nonzero_residual_improved_by_mean_maps=False,
    all_later_actual_source_families_and_finite_induction_verified=False,
    unconditional_all_pulse_family_bounds_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    actual_giant_profile_or_first_band_numerically_evaluated=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['independently_audit_all_later_source_families_and_finite_stage_induction_from_current_seed',
    'join_limiting_order_inequalities_stage_constants_and_full_potential_realization',
    'close_complete_physical_residual_viscosity_cutoffs_and_interactions_before_unconditional_claims']
DEPENDENCIES=('first_sourced_foundation','regular_source_foundation','physical_pulse_foundation','compact_mean_foundation',
    'real_signed_covariance_update','complex_signed_covariance_update','post_signed_mean_relaxation','complex_mean_relaxation')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='bd707aadc121c83563490495c412ab4c7d15e026',
        unchanged_parameters='Md64; unit viscosity; fixed h,N,positive lambda,t1=C_setup^-20; original radius roots supports carriers and operation order',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        current_input='B=rho=mu=1-kappa; C=3/2-kappa; D=H=C-kappa',
        generic_covariance_error='nu=min(rho+C-1/2,C+1/2-kappa,2C-1)=2-2*kappa; all three tie',
        mean_refresh='current signed state -> fast at H -> recompute pressure/defects -> slow at C -> recompute again',
        stale_target_obstruction='mu+H=5/2-3*kappa < C+1=5/2-kappa by 2*kappa',
        final_seed='nonzero D; full tangential mean nu=D+1/2; defects C+1; old error rho=B; Haar tangential 1 radial 2',
        pressure='all six W entries, full p_B including Pi1/Pi2, current p_m and whole derivative of Jz+c_rho*P',
        exact_tails='actual flat velocities retained in products; designated Gaussian tails additive, never inverse sources',
        scope='ONE current-input finite signed/mean join relative to exact core and finite lemmas; later families, full residual and blow-up remain open')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'first_sourced_foundation';sys.path.insert(0,str(p))
    try:return load('current_complex_foundation_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'first_sourced_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():history['first_sourced_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,history


def sources():
    paths=[HERE/n for n in ('README.md','complex_join_bounds.py','complex_join_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'first_sourced_foundation/sourced_checks.py',HERE.parent/'post_signed_mean_relaxation/relax.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(current_signed_and_mean_order_join=review(),**exact_reviews()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews();gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(actual_first_source_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Current complex foundation scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Later source/full residual obligations erased')
    keys=set(reviews())|{'actual_first_source_parent_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('A current complex foundation gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=925 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);args=p.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if args.verify_record:
        f=json.loads(args.verify_record.read_text());validate(f);compare_records(r,f)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('One current finite signed/mean join; later source families, full residual and blow-up remain open.')


if __name__=='__main__':main()
