#!/usr/bin/env python3
"""Regular source normalization, weighted edges and explicit tail-freezing correction."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json
from pathlib import Path
import sys
from source_bounds import acceptance_review,envelope_review,weighted_review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='regular-source-envelope-edge-estimates-with-explicit-freezing-correction'
DEPENDENCIES=('certified_tail_reconstruction','primary_covariance_foundation','sourced_inverse_foundation',
    'regular_source_foundation','compact_mean_foundation','canonical_finite_cycle_join')
FLAGS=dict(missing_displayed_direction_freezing_condition_identified=True,
    corrected_initial_tail_selector_contains_sixth_root_freezing_budget=True,
    length_budget_must_use_certified_frozen_margin=True,
    all_fixed_mixed_normalized_source_jets_have_Bell_majorants=True,
    generated_regular_source_grammar_has_full_path_and_entrance_inputs_relative_to_current_target_estimates=True,
    full_source_normal_moving_normal_and_Pi1_Pi2_pressure_retained=True,
    signed_reciprocal_root_half_weight_cancellation_quantified=True,
    two_edge_primitive_estimate_keeps_necessary_zero_total=True,
    actual_width_half_exponential_constants_and_physical_Haar_costs_retained=True,
    historical_sources_records_fields_and_forces_preserved=True,
    exact_parameter_and_jet_controls_separate_from_continuum_derivations=True,
    parent_displayed_tail_selector_alone_proves_direction_freezing=False,
    actual_construction_disproved_by_parameter_witness=False,
    historical_first_band_or_force_recertified=False,
    corrected_initial_choice_identified_as_old_realized_field_or_force=False,
    uniform_positive_unweighted_reciprocal_root_floor_at_shell_edges=False,
    every_actual_later_covariance_mean_target_verified_by_grammar_or_protocol=False,
    arbitrary_physically_flat_tail_has_forward_inverse_source_bound=False,
    actual_giant_fields_constants_or_first_band_numerically_evaluated=False,
    unconditional_pulse_complete_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['verify_actual_current_signed_target_full_zeta_jets_and_complete_covariance_mean_moment_identities_at_each_canonical_update_with_pressure_axial_velocity_and_reserved_hull',
    'use_corrected_freezing_budget_at_initial_tail_choice_before_reconstruction_not_as_an_identity_with_historical_fields_or_forces',
    'discharge_still_unverified_defining_map_and_complete_iteration_hypotheses_before_unconditional_claims']

def protocol():
    return dict(parent_commit='ffd0b394d7f0e5df08c9d171937ae3fb1d362d8e',
        correction='target freezing needs ell^6>=8*Cfreeze/gamma0, separate from Ls length bounds; historical evidence preserved',
        selector='Eplus=ceil(max(16,ell_old,512/h^2,sqrt(max(Lreq,Lcarrier)/(2r0)),(8*Cfreeze/gamma0)^(1/6)))',
        margin='Lreq evaluated at gamma0/2 or smaller; Cfreeze from actual smooth limiting direction not T0/zeta',
        normalization='P>0; nonnegative Bell recurrence bounds P*D^I(P^-1); quotient jets bounded by finite convolution with source jets',
        grammar='linear complete remainder; complete same-label curl products; mean times wave; finite harmonic n/m retained; additive flat tails excluded',
        edges='zeta=exp(-t1^2/ya^2-4/yb^2),t1=C_setup^-20; oriented primitive gains x^3/a near left and symmetric right only with zero total',
        signed='Xi has full zeta jets, reciprocal root loses zeta^(-1/2); delta_a<=C epsilon^(alpha-1) s0^(-1/2) sqrt(zeta) times finite edge powers',
        pressure='full n.K.t-nprime.t+n.f; actual Pi1/Pi2,axial means,all six covariance entries and moving pressure moment retained',
        remaining='actual current Xi weighted jets and both complete tangential moment identities with whole partial_Z(Jz+c_rho*Pmean)',
        scope='source functional estimates relative to current targets, corrected initial acceptance contract; no unconditional full cycle or identification with old force')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'certified_tail_reconstruction';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('source_edge_tail_parent',p/'audit.py');m=importlib.util.module_from_spec(s)
        sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old

@lru_cache(None)
def records():
    r={n:json.loads((HERE.parent/n/'evidence.json').read_text()) for n in DEPENDENCIES}
    parent_module().validate(r['certified_tail_reconstruction']);return r

@lru_cache(None)
def history():
    m=parent_module()
    for _ in range(3):m=m.parent_module()
    return m.provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','source_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    paths += [HERE.parent/'canonical_value_domain/acceptance.py',HERE.parent/'certified_tail_reconstruction/tail_bounds.py',HERE.parent/'regular_source_foundation/source_bounds.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def routing_review():
    r=records()
    return dict(checks=dict(
        parent_selector_did_not_display_Cfreeze='Cfreeze' not in r['certified_tail_reconstruction']['protocol']['domain'],
        old_first_band_still_not_verified=not r['certified_tail_reconstruction']['flags']['original_first_band_or_old_localized_force_acceptance_verified'],
        source_foundation_keeps_complete_pressure=r['regular_source_foundation']['flags']['complete_Pi1_Pi2_background_pressure_and_mean_moments_retained'],
        source_grammar_is_not_all_actual_target_verification=not r['regular_source_foundation']['flags']['all_actual_later_source_families_and_complex_targets_verified'],
        finite_cycle_requires_full_covariances_and_squares=r['canonical_finite_cycle_join']['flags']['complete_pressure_viscosity_all_six_covariances_and_nonlinear_squares_retained'],
        all_dependencies_remain_explicit=all(r[n]['remaining_obligations'] for n in DEPENDENCIES)))

@lru_cache(None)
def reviews():return dict(corrected_freezing_and_tail_contract=acceptance_review(),mixed_source_quotient_majorants=envelope_review(),
    reciprocal_roots_oriented_edges_and_moments=weighted_review(),current_operator_routing_and_limits=routing_review())
def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate();compare_records(expected,r)
    if not all(expected['gates'].values()):raise ValueError('Source/edge/freezing gate failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; historical evidence is not overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Functional source estimates and freezing correction only; actual current target/moment audit and unconditional proof remain.')
if __name__=='__main__':main()
