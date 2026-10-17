#!/usr/bin/env python3
"""All-order prescribed continuation maps through shape, relative to the exact core."""
import argparse
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from continuation_bounds import bounds_review,join_review
from continuation_checks import diagnostic_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='pre-modulation-prescribed-continuation-through-shape-all-fixed-jets'
FLAGS=dict(correct_integrated_reference_coordinates_not_flattened_physical_shear_used=True,
    every_fixed_jet_of_prescribed_continuation_through_shape_bounded_relative_to_exact_core=True,
    one_common_parameter_strip_for_every_fixed_radial_order=True,
    complete_reference_and_actual_five_moments_and_pressure_memory_retained=True,
    narrow_cutoff_derivative_costs_retained=True,
    finite_maps_have_no_unknown_parameter_derivative_feedback=True,
    historical_sources_evidence_and_parameters_unchanged=True,
    diagnostics_are_manufactured_not_actual_continuum_fields=True,
    exact_pressure_family_core_remains_an_explicit_input=True,
    first_collar_formulas_valid_after_native_reference_taper=False,
    pre_modulation_profile_equals_final_modulated_profile_on_whole_region=False,
    uniform_in_C_or_all_derivative_orders_constants_verified=False,
    all_order_finite_N_modulation_phase_and_exact_repair_independently_verified=False,
    global_backward_pressure_identity_reproved_here=False,
    unconditional_actual_pulse_bounds_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['audit_same_finite_N_modulation_phase_inverse_and_exact_five_row_repair_at_all_fixed_orders',
    'join_annular_and_exterior_profile_maps_complete_global_pressure_moments_weighted_edges_and_reserved_supports',
    'only_then_promote_background_actual_pulse_mean_and_complete_residual_bounds_beyond_conditional_claims']
DEPENDENCIES=('activated_collar_foundation','core_all_order_foundation','reduced_continuation_margins',
    'shape_join_pressure','compact_jet_envelope','repair_state_bounds')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='bde6e3d9d212a706b54d8688941379400f60f79d',
        fixed_width='t1=kappa0=C^-20',final_eta_strip='rho/(12800*S)',
        domain='0<=X<=Xsep=110*exp(2000*B), pre-modulation prescribed profile',
        reference_coordinates='original five forward moments and full Pi_axis+Cp, not flattened physical shear',
        actual_maps='known reference p1_r,ns_r -> prescribed slopes -> primitives and exponentials -> actual five moments -> shape',
        no_parameter_support_or_frozen_evidence_change=True,
        numerical_diagnostic='independent exact rational bivariate jets plus legacy raw second jets',
        scope='relative finite all-order seminorm theorem; final modulation/repair and global foundation remain open')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'activated_collar_foundation';sys.path.insert(0,str(p))
    try:return load('prescribed_continuation_activated_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'activated_collar_foundation/evidence.json').read_text())
    m.validate(f);m.compare_records(r,f)
    inherited=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():
        inherited['activated_collar_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,inherited


def sources():
    paths=[HERE/n for n in ('README.md','continuation_bounds.py','continuation_checks.py','audit.py',
        'test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    # Also bind the two independent executable comparison sources directly.
    paths += [HERE.parent/'core_all_order_foundation/coordinate_checks.py',HERE.parent/'compact_jet_envelope/jets.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return dict(finite_seminorm_and_parameter_domain=bounds_review(),
    defining_fields_full_memory_and_smooth_shape_join=join_review(),
    independent_rational_prescribed_flow_and_raw_pressure_controls=diagnostic_review())


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    parent,inherited=parent_review();r=reviews()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(exact_core_actual_collar_and_conditional_chain_reproduced=all(parent['gates'].values()),
        historical_source_bytes_preserved=bool(inherited) and all(inherited.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():raise ValueError('Continuation scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Global foundation gate erased')
    expected=reviews();keys=set(expected)|{'exact_core_actual_collar_and_conditional_chain_reproduced','historical_source_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('A continuation gate failed')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=785 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
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
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} source/record bindings and {len(r['inherited_provenance'])} inherited checks pass")
    print('PRE-MODULATION continuation relative to exact core only; final global profile, pulses, residual and blow-up remain open.')


if __name__=='__main__':main()
