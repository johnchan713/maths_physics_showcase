#!/usr/bin/env python3
"""Join current-source foundations to every finite complete correction cycle."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from finite_bounds import seed_review,half_line_review,physical_review
from finite_checks import reviews as exact_reviews

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-source-complete-finite-induction-relative-to-core-and-finite-analytic-lemmas'
FLAGS=dict(current_complex_seed_B0_three_halves_minus_two_kappa_joined=True,
    all_finite_source_expressions_closed_by_supported_complete_curl_grammar=True,
    complete_finite_cycle_gain_one_half_minus_kappa_on_whole_half_line=True,
    seven_nonzero_and_two_covariance_rows_limit_available_gain=True,
    complete_pressure_viscosity_cutoffs_and_all_interactions_retained=True,
    full_updated_covariance_and_mean_pressure_recomputed_after_each_change=True,
    fresh_postfast_defects_actual_A0_A1_and_old_Haar_means_retained=True,
    both_background_pressure_tails_and_whole_moving_pressure_product_retained=True,
    all_fixed_jets_all_accepted_bands_and_smooth_endpoint_supports_joined=True,
    finite_harmonics_and_stage_constants_may_grow_without_retuning_domain=True,
    exact_additive_Gaussian_tails_never_forward_inverted=True,
    each_fixed_target_power_and_jet_order_has_some_finite_state=True,
    retained_exact_core_and_finite_defining_inverse_gap_lemmas_remain_dependencies=True,
    historical_evidence_geometry_parameters_roots_width_and_operation_order_unchanged=True,
    manufactured_finite_mode_checks_not_continuum_suprema=True,
    stage_constants_uniform_or_fixed_q_contraction_verified=False,
    arbitrary_stage_j_of_q_substitution_justified=False,
    one_infinite_realized_field_verified_from_current_foundations=False,
    downstream_localization_force_and_comparison_promoted_unconditionally=False,
    all_retained_analytic_dependencies_discharged=False,
    unconditional_complete_residual_improvement_verified=False,
    actual_giant_or_thousand_stage_fields_numerically_evaluated=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['independently_join_complete_potentials_direct_swirl_pressure_gauges_and_diagonal_cutoffs_to_one_current_realization',
    'then_join_endpoint_localization_admissible_force_energy_and_global_comparison_to_that_same_field',
    'discharge_retained_exact_core_and_finite_analytic_existence_inverse_gap_dependencies_before_unconditional_claims']
DEPENDENCIES=('current_complex_foundation','first_sourced_foundation','regular_source_foundation','physical_pulse_foundation',
    'compact_mean_foundation','finite_stage_induction','current_potential_realization','residual_cycle_closure')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='89f9c856d8800a71f16c36e23819babb61aa0c63',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        unchanged_parameters='Md64; unit viscosity; original h,N,positive lambda,width C_setup^-20,roots,geometry,carriers and operation order',
        current_seed='B0=3/2-2*kappa; full mean B0+1/2; defects 5/2-kappa; rho=mu=1-kappa; Haar tangent 1 radial 2',
        all_finite_orders='B_j=3/2-2*kappa+j*(1/2-kappa), C_j=B_j+1/2',
        gain_obstruction='s<=min(1/2-kappa,rho-1/2,mu-1/2)=1/2-kappa',
        proof_domain='entire B>=B0 by affine slopes and boundary margins, supported-source expression induction and retained analytic maps',
        cycle='regular sourced wave; current generic signed covariance; fresh fast mean; current slow five-map',
        pressure='all six covariance entries, axial mean velocity and pressure, both Pi1/Pi2, whole derivative of Jz+c_rho*P',
        source_split='regular supported expression only; designated Gaussian linear tails additive never forward-inverted; actual flat velocities stay in products',
        harmonics='M_j<=4*2^j; finite stage constants may grow; same qstar and enlarged charts 4Q<qstar',
        physical_bound='for each fixed finite j,m: |R(U_j)|_m<=C_jm*q^(h*B_j-(2*A+3/2+2*m))',
        quantifier='for each fixed target and jet there exists a finite state; one simultaneous realized field requires a separate diagonal argument',
        scope='all-finite current-source foundation join relative to retained exact-core and finite analytic lemmas; no unconditional realization force or blow-up')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'current_complex_foundation';sys.path.insert(0,str(p))
    try:return load('finite_foundation_current_complex_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'current_complex_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():history['current_complex_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,history


def sources():
    paths=[HERE/n for n in ('README.md','finite_bounds.py','finite_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'regular_source_foundation/source_bounds.py',HERE.parent/'current_complex_foundation/complex_join_bounds.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(current_actual_seed=seed_review(parent_review()[0]),whole_half_line_cycle=half_line_review(),
    complete_physical_residual_and_quantifiers=physical_review(),**exact_reviews()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews();gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(current_complex_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Finite-source foundation scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Realization/analytic dependency obligations erased')
    keys=set(reviews())|{'current_complex_parent_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('A finite-source foundation gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=943 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);args=p.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if args.verify_record:
        f=json.loads(args.verify_record.read_text());validate(f);compare_records(r,f)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Current all-finite source join; one-field realization, retained analytic dependencies and blow-up remain open.')


if __name__=='__main__':main()
