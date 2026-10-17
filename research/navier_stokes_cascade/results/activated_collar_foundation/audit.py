#!/usr/bin/env python3
"""Actual flat collar common-strip certificate relative to the exact pressure-family core."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from activation_bounds import domain_review,formula_review,width_review,collar_review
from activation_checks import diagnostic_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='actual-flat-activation-common-parameter-collar-and-all-order-leading-jets'
FLAGS=dict(actual_activated_fields_not_unmodified_seed_are_input=True,
    common_eta_strip_for_all_fixed_radial_orders_derived=True,
    full_five_moment_and_pressure_memory_preserved=True,
    first_activated_background_collar_certificate_supplied_relative_to_exact_core=True,
    all_fixed_physical_scale_leading_jet_bounds_extended_to_actual_first_collar=True,
    fixed_width_derivative_costs_and_difference_gains_explicit=True,
    exact_pressure_family_core_and_inherited_parameters_remain_explicit_hypotheses=True,
    numerical_diagnostics_not_actual_enormous_fields_or_continuum_certificates=True,
    actual_profile_radially_analytic_at_flat_join=False,
    radial_constants_uniform_as_width_or_C_vary=False,
    later_global_profile_modulation_and_repair_independently_verified=False,
    completed_global_backward_pressure_identity_reproved_here=False,
    unconditional_background_or_actual_pulse_bounds_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['independently_audit_later_reference_continuation_shape_and_compact_modulation_repair_on_the_same_actual_profile',
    'independently_verify_completed_global_pressure_moments_weighted_edges_and_reserved_patch_inputs',
    'then_join_the_common_collar_to_all_order_background_inverse_and_actual_pulse_bounds_before_unconditional_claims']
DEPENDENCIES=('core_all_order_foundation','inner_collar_activation','axis_core_attachment',
    'core_growth_activation','background_jet_foundation','coefficient_support_audit')


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='b1c66663d020a3920aafc58cb46aae421aeda8c1',
        input='exact infinite pressure-family core on Y<=4.1; actual final leading field unchanged only on Y<=4',
        fixed_width='t1=kappa0=C^-20',collar='0<=y=log(X/Xa)<=tc/4, tc=t1/(10*sqrt(log C)), Xa=4/Lambda',
        final_eta_strip='rho/(800*S)',parameter_holomorphy_not_radial_analyticity=True,
        field_difference_gain='t1^(1-r)',moment_and_pressure_difference_gain='t1^(2-r)',
        coarse_full_radial_cost='t1^(-r)=C^(20*r); constants depend on frozen seed and order',
        actual_pressure='same Pi_axis and core trace plus actual X*F_a^2 forward increment',
        no_parameter_or_support_change=True,diagnostic_seed='historical manufactured a=3 analytic seed, not actual core',
        diagnostic_widths=['.01','.001','.0001'],diagnostic_s='.24',diagnostic_eta='.3',diagnostic_kappa0='.1',
        diagnostic_radial_order=4,diagnostic_eta_orders=list(range(7)),diagnostic_digits=75,
        quadrature_orders=[32,48],historical_formula_order=24,
        scope='common actual first collar relative to exact core; later global profile, backward pressure, full background/pulses, residual and blow-up not independently verified')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'core_all_order_foundation';sys.path.insert(0,str(p))
    try:return load('activated_collar_core_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'core_all_order_foundation/evidence.json').read_text())
    m.validate(f);m.compare_records(r,f)
    inherited=dict(r['inherited_provenance'])
    for path,digest in f['source_hashes'].items():
        inherited['core_all_order_foundation:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return r,inherited


def sources():
    names=('README.md','activation_bounds.py','activation_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')
    paths=[HERE/n for n in names]+[HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():
    return dict(actual_domain_and_cutoff_order=domain_review(),complete_actual_forward_formulas=formula_review(),
        all_order_width_losses_and_log_conversion=width_review(),common_collar_and_physical_transfer=collar_review(),
        higher_jet_complex_parameter_and_full_memory_controls=diagnostic_review())


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    parent,inherited=parent_review();r=reviews()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates['infinite_core_and_current_conditional_chain_reproduced']=all(parent['gates'].values())
    gates['historical_source_bytes_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():raise ValueError('Actual collar scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Global foundation gate erased')
    expected=reviews();keys=set(expected)|{'infinite_core_and_current_conditional_chain_reproduced','historical_source_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('An activated collar gate failed')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=771 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
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
    print(f"{STATUS}: {len(r['gates'])} gates and {len(r['inherited_provenance'])} inherited hashes pass")
    print('Actual first collar relative to exact pressure-family core only; global foundation and unconditional blow-up remain open.')


if __name__=='__main__':main()
