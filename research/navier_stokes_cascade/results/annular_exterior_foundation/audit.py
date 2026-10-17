#!/usr/bin/env python3
"""Actual annular/exterior smooth-input and all-order global moment foundation."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from exterior_bounds import foundation_review,tail_review,support_review
from exterior_checks import diagnostic_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='annular-exterior-all-order-leading-regularity-and-global-moment-foundation-relative-to-exact-roots'
FLAGS=dict(whole_compact_smooth_input_to_same_finite_N_maps_supplied_relative_to_exact_core_and_finite_existence_bounds=True,
    eta_dependent_annular_Q_derivatives_included_at_every_fixed_order=True,
    same_accepted_annular_and_exterior_roots_have_all_fixed_jets=True,
    Gamma_heat_endpoint_and_integrable_global_moment_majorants_all_fixed_orders_derived=True,
    same_exact_five_global_moment_functions_and_all_eta_derivatives_identified=True,
    actual_reserved_M_equals_mu_eta_and_V_equals_minus_mu_not_zero=True,
    complete_intermediate_pressure_memory_and_exceptional_pressure_orders_retained=True,
    common_analytic_FIRST_collar_preserved_and_final_leading_closed_rectangle_regularities_joined=True,
    all_fixed_physical_leading_jet_bounds_on_finite_profile_ranges_joined=True,
    historical_sources_evidence_parameters_and_supports_preserved=True,
    diagnostics_are_manufactured_not_actual_continuum_suprema=True,
    finite_order_existence_inverse_and_gap_lemmas_remain_explicit_dependencies=True,
    heat_tail_analytic_in_eta_at_poles=False,
    estimates_uniform_in_h_to_zero_or_all_orders_or_N=False,
    all_order_weighted_T1_over_T0_and_background_edge_join_independently_verified=False,
    unconditional_actual_pulse_bounds_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['independently_audit_exact_flat_edge_leading_and_first_stress_quotient_and_weighted_background_join',
    'join_formal_coefficient_inverse_and_supports_to_actual_leading_inputs_retaining_I4_pressure_tails',
    'verify_uniform_actual_physical_pulses_complete_covariance_means_and_complete_residual_before_unconditional_claims']
DEPENDENCIES=('modulation_repair_all_order','axis_matching_audit','shape_join_pressure',
    'pulse_moment_audit','heat_exterior_audit','global_moment_stress','coefficient_support_audit','regional_cone_join')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='760ee0033020260296c1f6f3f3a264cfc0f4402c',
        unchanged_frequency='N=1+floor(H^32)',annular_map='B constant, Q_eta U-square coefficient exp(6/5)*(1+eta^2)^2/P^2',
        heat_domain='Z=2*(1-eta^2)/X>=0, fixed h>0, one-sided pole jets not global analyticity',
        reserved_fields='U0=0,M0=mu*eta,V0=-mu; e3 positive; full K3(eta) pressure memory retained',
        input_join='whole compact smooth pre-modulation input derived from exact core and finite defining maps',
        no_parameter_support_or_historical_evidence_change=True,
        scope='all-order leading regularity/global moment foundation relative to exact core and retained finite existence/gap bounds; weighted first-stress and physical pulse/residual foundation open')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'modulation_repair_all_order';sys.path.insert(0,str(p))
    try:return load('annular_exterior_all_order_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'modulation_repair_all_order/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    inherited=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():
        inherited['modulation_repair_all_order:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,inherited


def sources():
    paths=[HERE/n for n in ('README.md','exterior_bounds.py','exterior_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_all_order_foundation/coordinate_checks.py',HERE.parent/'modulation_repair_all_order/map_checks.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(nonconstant_annular_and_finite_exterior_maps=foundation_review(),
    all_order_heat_and_global_moment_domination=tail_review(),actual_reserved_plateau_and_pressure_memory=support_review(),
    independent_nonconstant_root_and_heat_controls=diagnostic_review()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    parent,inherited=parent_review();r=reviews()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(parent_map_foundation_and_retained_conditional_chain_reproduced=all(parent['gates'].values()),
        historical_source_bytes_preserved=bool(inherited) and all(inherited.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():raise ValueError('Annular/exterior foundation scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Weighted-edge/physical gate erased')
    expected=reviews();keys=set(expected)|{'parent_map_foundation_and_retained_conditional_chain_reproduced','historical_source_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('An exterior foundation gate failed')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=818 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(expected,r.get('reviews',{}),'root.reviews')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    args=p.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; historical frozen evidence is never overwritten')
    r=generate();validate(r)
    if args.verify_record:
        f=json.loads(args.verify_record.read_text());validate(f);compare_records(r,f)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings and {len(r['inherited_provenance'])} inherited checks pass")
    print('Retained exact core and finite existence lemmas remain inputs; weighted stress, actual pulses, full residual and blow-up are not promoted.')


if __name__=='__main__':main()
