#!/usr/bin/env python3
"""Actual weighted first-stress/edge and background-input foundation join."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from edge_bounds import actual_input_review,outer_review,background_review
from edge_checks import diagnostic_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='actual-weighted-first-stress-and-background-foundation-joined-relative-to-exact-core-and-finite-lemmas'
FLAGS=dict(actual_leading_inputs_to_formal_coefficient_maps_supplied=True,
    all_fixed_weighted_leading_and_first_stress_edges_independently_joined=True,
    T1_over_abs_T0_smooth_on_closed_annulus_including_eta_zero=True,
    actual_Gamma_eta_zero_cubic_quotient_Xb_over_eight_derived=True,
    exact_common_higher_stress_supports_and_backward_moments_retained=True,
    one_admissible_diagonal_schedule_controls_all_fixed_relative_background_jets=True,
    growth_and_oriented_background_cone_on_one_small_q_domain_derived=True,
    full_pressure_memory_I4_first_two_orders_and_radial_velocity_retained=True,
    exact_pressure_family_core_and_finite_existence_inverse_gap_lemmas_remain_dependencies=True,
    historical_evidence_sources_parameters_width_frequency_and_supports_preserved=True,
    finite_diagnostics_not_continuum_suprema=True,
    heat_pole_Taylor_series_assumed_convergent=False,
    constants_uniform_in_h_order_N_or_actual_coefficient_index=False,
    actual_giant_fields_or_diagonal_scales_numerically_evaluated=False,
    unconditional_uniform_actual_physical_pulses_verified=False,
    unconditional_complete_covariance_means_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['independently_join_actual_uniform_physical_pulse_values_derivatives_and_supports_at_every_scale',
    'audit_complete_wave_covariance_means_with_axial_velocity_radial_flow_and_Pi1_Pi2',
    'verify_complete_residual_with_viscosity_all_cutoffs_and_interactions_before_unconditional_claims']
DEPENDENCIES=('annular_exterior_foundation','weighted_background_edges','inner_collar_activation',
    'outer_edge_stability','coefficient_support_audit','background_jet_foundation','global_moment_stress','regional_cone_join')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='a779312a3036c03ca7abea540b18e19785c9bc03',
        unchanged_parameters='Md64; t1=C^-20; fixed h and N=1+floor(H^32)',
        exact_weight='exp(-t1^2/ya^2-4/yb^2)',total_order_loss='3*m+3',
        first_stress_quotient='(b1/sqrt(b0^2+yb^12*bz^2),0)',
        eta_zero_cubic='W1theta(delta,0)=Xb*delta^3/8+O(delta^4)',
        derivative_family='Euler q*d_q at fixed X,eta; partial_X; partial_eta',
        one_diagonal_schedule_all_earlier_constraints_and_new_finite_quotient_seminorms=True,
        retain_exact_core_and_finite_lemmas_not_manuscript_conclusion=True,
        scope='weighted leading/first-stress and normalized background foundation; actual pulses complete means residual and blow-up not promoted')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'annular_exterior_foundation';sys.path.insert(0,str(p))
    try:return load('weighted_stress_foundation_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'annular_exterior_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    inherited=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():
        inherited['annular_exterior_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,inherited


def sources():
    paths=[HERE/n for n in ('README.md','edge_bounds.py','edge_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'outer_edge_stability/tail.py',HERE.parent/'weighted_background_edges/edges.py',HERE.parent/'background_jet_foundation/foundation.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(actual_inner_factor_and_coefficient_support_inputs=actual_input_review(),
    all_order_outer_factor_and_eta_zero_cubic=outer_review(),
    single_schedule_weighted_background_and_common_cone=background_review(),
    independent_exact_kernel_heat_and_budget_controls=diagnostic_review()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,inherited=parent_review();r=reviews()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(actual_leading_parent_and_retained_conditional_chain_reproduced=all(p['gates'].values()),
        historical_bytes_preserved=bool(inherited) and all(inherited.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():raise ValueError('Weighted foundation scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Actual pulse/mean/residual gates erased')
    expected=reviews();keys=set(expected)|{'actual_leading_parent_and_retained_conditional_chain_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('A weighted foundation gate failed')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=836 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(expected,r.get('reviews',{}),'root.reviews')
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
    print('Exact core and retained finite lemmas remain inputs; actual pulses, complete means/residual and blow-up remain open.')


if __name__=='__main__':main()
