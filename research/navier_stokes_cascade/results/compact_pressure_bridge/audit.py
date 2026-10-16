#!/usr/bin/env python3
"""Replay the whole-activation and conditional compact pulse bridge."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from bridge import activation_review, memory_review, operator_review, pulse_review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'whole-activation-and-conditional-compact-pulse-bridge-audited'
FLAGS = dict(whole_first_activation_relaxed_cone_bound_derived=True,
             exact_equal_field_pressure_memory_identity_derived=True,
             exact_five_moment_restoration_implies_forward_outer_equality=True,
             quantitative_compact_pulse_values_conditional_on_input_margins=True,
             post_activation_continuation_margins_independently_verified=False,
             completed_global_backward_pressure_normalization_verified=False,
             whole_leading_profile_independently_verified=False,
             unconditional_uniform_actual_pulse_estimate_verified=False,
             unconditional_complete_residual_improvement_verified=False,
             blowup_verified=False)
REMAINING = ['audit_actual_post_activation_reduced_shear_continuation_and_both_final_tapers',
             'audit_shape_join_outer_dependencies_and_completed_backward_pressure_normalization',
             'verify_coefficient_and_reserved_mean_supports',
             'verify_all_scale_actual_pulse_derivatives_and_supports',
             'retain_axial_velocity_and_pressure_in_complete_covariance_mean_corrections',
             'prove_complete_residual_gain_with_viscosity_cutoffs_and_all_interactions']


def protocol():
    return dict(parent_commit='08c4c30d8a9f8de4ebd00b7484db6ab55530fe24',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                same_fixed_parameters_and_axis_datum=True,
                activation_s_interval=['0','1'], unchanged_width='C^-20',
                activation_Pc_lower='11/5', activation_squared_cone_loss_upper='C^-30',
                original_compact_input_margins_remain_hypotheses=True,
                compact_R='4*G^3', compact_kappa='G^-4', compact_Fmin='C^-11',
                compact_squared_cone_deficit='kappa/(64*R^6)',
                pressure_memory_rows=['M','I','J','S','Cp'],
                local_forward_pressure_equality_is_not_global_backward_normalization=True,
                diagnostic_operator_digits=50, diagnostic_operator_orders=[48,72],
                fixtures_are_not_actual_core_or_actual_pulses=True,
                diagnostics_are_not_continuum_certificates=True)


def inherited_hashes():
    parent = HERE.parent/'core_growth_activation'
    spec = importlib.util.spec_from_file_location('compact_bridge_frozen_parent',parent/'audit.py')
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(parent))
    try:
        spec.loader.exec_module(module)
        checks = module.inherited_hashes()
    finally:
        sys.path.pop(0)
    record = json.loads((parent/'evidence.json').read_text())
    for path, expected in record['source_hashes'].items():
        checks['core_growth_activation:'+path] = hashlib.sha256((PROJECT/path).read_bytes()).hexdigest() == expected
    return checks


def source_hashes():
    files = [HERE/n for n in ('README.md','bridge.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/'evidence.json' for n in ('core_growth_activation','loop_modulation_bounds','repair_state_bounds','outer_edge_stability')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != protocol():
        raise ValueError('Protocol changed')
    reviews = dict(whole_activation_scalar_budget=activation_review(),
                   exact_pressure_memory=memory_review(),
                   compact_pulse_value_identities=pulse_review(),
                   positive_whole_interval_operator=operator_review())
    inherited = inherited_hashes()
    gates = {k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews.items()}
    gates['frozen_provenance_preserved'] = bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
                inherited_provenance=inherited,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def validate(record):
    if record.get('status') != STATUS or record.get('protocol') != protocol() or record.get('flags') != FLAGS:
        raise ValueError('Changed derivation scope or parameters')
    if record.get('remaining_obligations') != REMAINING:
        raise ValueError('An outstanding obligation was removed')
    names = {'whole_activation_scalar_budget','exact_pressure_memory',
             'compact_pulse_value_identities','positive_whole_interval_operator'}
    if set(record.get('reviews',{})) != names:
        raise ValueError('A review is missing')
    if set(record.get('gates',{})) != names|{'frozen_provenance_preserved'} or not all(record['gates'].values()):
        raise ValueError('A bridge gate failed')
    if not record.get('inherited_provenance') or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited bytes changed')
    if not all(v.get('checks') and all(v['checks'].values()) for v in record['reviews'].values()):
        raise ValueError('A diagnostic failed')


def compare_records(a,b,path='root'):
    if type(a) is not type(b):
        raise ValueError('Changed type at '+path)
    if isinstance(a,dict):
        if set(a) != set(b):
            raise ValueError('Changed keys at '+path)
        for k in a:
            compare_records(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if len(a) != len(b):
            raise ValueError('Changed list at '+path)
        for i,(x,y) in enumerate(zip(a,b)):
            compare_records(x,y,path+'.'+str(i))
    elif isinstance(a,float):
        if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-40):
            raise ValueError('Changed diagnostic at '+path)
    elif a != b:
        raise ValueError('Changed exact value at '+path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit('Output exists; choose a fresh path')
    result = generate()
    validate(result)
    if args.verify_record:
        frozen = json.loads(args.verify_record.read_text())
        validate(frozen)
        compare_records(result,frozen)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(STATUS)
    print(f"{len(result['gates'])} gates pass; {len(result['inherited_provenance'])} inherited hashes agree.")
    print('Whole activation derived; compact pulse values conditional; global pressure and full residual open.')


if __name__ == '__main__':
    main()
