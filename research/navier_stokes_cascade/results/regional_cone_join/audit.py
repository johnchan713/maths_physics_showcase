#!/usr/bin/env python3
"""Replay corrected regional input gaps and the leading value/cone join."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import mpmath as mp
from cones import comparison_review, equation_review, joined_review, outward_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='regional-input-gaps-and-leading-value-cone-joined'
FLAGS=dict(original_compact_regional_value_margins_derived=True,
           first_outer_ramp_signed_comparison_corrected=True,
           compact_input_cone_hypotheses_discharged=True,
           same_finite_frequency_leading_value_cone_joined=True,
           exact_global_moment_ledger_applies_to_leading_value_construction=True,
           relies_on_replayed_analytic_core_envelope_loop_and_inverse_lemmas=True,
           all_order_leading_profile_certificate_independently_verified=False,
           positive_order_and_reserved_mean_supports_verified=False,
           unconditional_uniform_actual_pulse_estimate_verified=False,
           unconditional_complete_residual_improvement_verified=False,
           blowup_verified=False)
REMAINING=['verify_positive_order_coefficient_support_and_untouched_reserved_mean_intervals',
           'audit_all_order_weighted_edges_required_by_the_background_jet_certificate',
           'verify_all_scale_actual_pulse_values_derivatives_and_supports',
           'retain_axial_velocity_and_pressure_in_complete_covariance_mean_corrections',
           'prove_complete_residual_gain_with_viscosity_cutoffs_and_all_interactions']
DEPENDENCIES=('global_moment_stress','compact_jet_envelope','loop_modulation_bounds',
              'repair_state_bounds','compact_pressure_bridge','reduced_continuation_margins',
              'shape_join_pressure','axial_stress_audit','intermediate_decay_audit',
              'heat_exterior_audit','outer_edge_stability')


def protocol():
    return dict(parent_commit='1cbfe7e7ca29321ce510ac6e07e5e24be9b9015c',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                same_parameters_widths_roots_supports_and_frequency=True,
                selected_family='Md64, T=exp(64)+10, lambda=exp(-4T), h=lambda^2',
                original_required_gap_floor='C^-1024 > Aenv^-1',Aenv='C^4096',
                original_a_floor='C^-22',original_Pc_minus_two_floor='.1',
                original_third_fourth_required_only_where_vs_at_least_two=True,
                corrected_first_ramp_floor='1.1*exp(-.6)-(.01/.98)*(exp(1)-1) > .586 > .5',
                old_outer_ramp_stronger_kernel_inference_superseded=True,
                negative_control_is_comparison_model_not_actual_ramp=True,
                left_collar_vs_floor='2.29',right_collar_vs_minus_two='2*lambda > C^-1',
                compact_final_gap_floor='kappa/2, kappa=G^-4',
                accepted_N='1+floor(Hhier^32)',
                compact_right_equals_heat_compensation_left=True,
                later_cone_uses_heat_compensated_profile_with_exact_pressure_memory=True,
                diagnostic_interval_digits=[100,160],
                actual_frequency_or_profile_numerically_resolved=False,
                leading_value_cone_is_not_all_scale_physical_pulse_or_full_residual=True)


def parent_review():
    parent=HERE.parent/'global_moment_stress'
    spec=importlib.util.spec_from_file_location('regional_frozen_global_audit',parent/'audit.py')
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(parent))
    try:
        spec.loader.exec_module(module)
        generated=module.generate();module.validate(generated)
        frozen=json.loads((parent/'evidence.json').read_text())
        module.compare_records(generated,frozen)
    finally:
        sys.path.pop(0)
    inherited=dict(generated['inherited_provenance'])
    for path,expected in frozen['source_hashes'].items():
        inherited['global_moment_stress:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==expected
    return generated,inherited


def source_hashes():
    paths=[HERE/n for n in ('README.md','cones.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():
        raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review()
    reviews=dict(exact_finite_cone_identities=equation_review(),
                 corrected_signed_outer_ramp=comparison_review(),
                 actual_regional_input_and_acceptance_budgets=outward_review(),
                 compact_and_compensated_exterior_join=joined_review())
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews.items()}
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    gates['exact_global_moment_dependency_replayed']=all(parent['gates'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
                inherited_provenance=inherited,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def validate(record):
    if record.get('status')!=STATUS or record.get('protocol')!=protocol() or record.get('flags')!=FLAGS:
        raise ValueError('Changed fixed choices or scientific scope')
    if record.get('remaining_obligations')!=REMAINING:
        raise ValueError('An outstanding obligation was removed')
    names={'exact_finite_cone_identities','corrected_signed_outer_ramp',
           'actual_regional_input_and_acceptance_budgets','compact_and_compensated_exterior_join'}
    if set(record.get('reviews',{}))!=names or set(record.get('gates',{}))!=names|{'frozen_provenance_preserved','exact_global_moment_dependency_replayed'}:
        raise ValueError('A review or gate was removed')
    if not all(record['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in record['reviews'].values()):
        raise ValueError('A regional cone/join check failed')
    if not record.get('inherited_provenance') or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited bytes changed')
    actual=record['reviews']['compact_and_compensated_exterior_join']
    expected=joined_review()
    if any(actual.get(k)!=expected[k] for k in ('regions','input_floor','compact_final_four_gap_floor','selected_N','immediate_exterior')):
        raise ValueError('Changed endpoint coverage or finite-frequency join')


def compare_records(a,b,path='root'):
    if type(a) is not type(b):
        raise ValueError('Changed type at '+path)
    if isinstance(a,dict):
        if set(a)!=set(b):
            raise ValueError('Changed keys at '+path)
        for k in a:
            compare_records(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if len(a)!=len(b):
            raise ValueError('Changed list at '+path)
        for i,(x,y) in enumerate(zip(a,b)):
            compare_records(x,y,path+'.'+str(i))
    elif isinstance(a,float):
        if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-70):
            raise ValueError('Changed numerical diagnostic at '+path)
    elif isinstance(a,str) and a!=b and path.rsplit('.',1)[-1] in ('lower','upper'):
        with mp.workdps(170):
            x,y=mp.mpf(a),mp.mpf(b)
            if not mp.isfinite(x) or not mp.isfinite(y) or abs(x-y)>mp.mpf('1e-55')*max(1,abs(x),abs(y)):
                raise ValueError('Changed interval endpoint at '+path)
    elif a!=b:
        raise ValueError('Changed exact value at '+path)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit('Output exists; choose a fresh path')
    record=generate();validate(record)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text());validate(frozen);compare_records(record,frozen)
    args.output.write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
    print(STATUS)
    print(f"{len(record['gates'])} gates pass; {len(record['inherited_provenance'])} inherited hashes agree.")
    print('Leading value/cone joined using analytic lemmas; all-order edges, supports, physical pulses and full residual remain open.')


if __name__=='__main__':
    main()
