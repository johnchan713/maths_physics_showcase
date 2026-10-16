#!/usr/bin/env python3
"""Replay global leading moments, signed stresses and exact exterior defects."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from ledger import diagnostic_review, exact_review, ledger_review, scaling_review, tail_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='exact-leading-global-moments-and-signed-stress-identification-audited'
FLAGS=dict(global_leading_moments_for_exact_completed_roots_derived=True,
           exact_five_functional_exterior_defect_map_derived=True,
           signed_weighted_leading_residual_totals_zero_for_exact_roots=True,
           outer_forward_backward_leading_stress_identification_for_exact_roots=True,
           compact_input_hypotheses_remain_dependencies=True,
           remaining_outer_regional_cone_inputs_independently_verified=False,
           whole_leading_profile_independently_verified=False,
           coefficient_and_reserved_mean_supports_verified=False,
           unconditional_uniform_actual_pulse_estimate_verified=False,
           unconditional_complete_residual_improvement_verified=False,
           blowup_verified=False)
REMAINING=['audit_outer_regional_cone_inputs_jointly_with_inner_shape_join_margins',
           'verify_positive_order_coefficient_and_reserved_mean_supports',
           'verify_all_scale_actual_pulse_values_derivatives_and_supports',
           'retain_axial_velocity_and_pressure_in_complete_covariance_mean_corrections',
           'prove_complete_residual_gain_with_viscosity_cutoffs_and_all_interactions']


def protocol():
    return dict(parent_commit='cf7a30cb7dddd82278f9e47031ad97458e3eed45',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                same_fixed_parameters_axis_datum_roots_and_supports=True,
                source_equations=['4.7','4.11','4.15-4.16','A.11','A.34-A.46'],
                global_rows=['M','renormalized_I','J','S','Cp'],
                global_targets=['0','0','0','0','-Pi0'],
                renormalized_I='integral_0^infinity(H-Hpow), Hpow=k*X^-h, k independent of eta',
                pressure_defect='Pi0+Cp(infinity)',
                pure_heat_only_after_all_leading_edits=True,
                completed_exact_root_and_compact_input_hypotheses_explicit=True,
                leading_stress_not_complete_PDE_residual=True,
                positive_order_support_premise_not_closed=True,
                diagnostic_h='1/8, not selected actual h',diagnostic_digits=[50,70],
                diagnostic_quadrature_orders=[24,36],diagnostic_Taylor_degree=12,
                diagnostic_small_Z_threshold='1e-6',
                diagnostics_are_not_actual_profile_supremum_or_exact_root_proofs=True)


def parent_review():
    parent=HERE.parent/'shape_join_pressure'
    spec=importlib.util.spec_from_file_location('global_moment_frozen_shape_audit',parent/'audit.py')
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(parent))
    try:
        spec.loader.exec_module(module)
        generated=module.generate();module.validate(generated)
        frozen=json.loads((parent/'evidence.json').read_text())
        module.compare_records(generated,frozen)
        inherited=module.inherited_hashes()
    finally:
        sys.path.pop(0)
    for path,expected in frozen['source_hashes'].items():
        inherited['shape_join_pressure:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==expected
    return generated,inherited


def source_hashes():
    files=[HERE/n for n in ('README.md','ledger.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/'evidence.json' for n in ('shape_join_pressure','pulse_moment_audit','heat_exterior_audit','repair_state_bounds','outer_edge_stability')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():
        raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review()
    reviews=dict(exact_exterior_defects=exact_review(),signed_edit_ledger=ledger_review(),
                 integrable_renormalized_tail=tail_review(),physical_signed_scaling=scaling_review(),
                 labelled_heat_and_Fubini_diagnostics=diagnostic_review())
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews.items()}
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    gates['exact_root_dependency_acceptance']=all(parent['gates'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
                inherited_provenance=inherited,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def validate(record):
    if record.get('status')!=STATUS or record.get('protocol')!=protocol() or record.get('flags')!=FLAGS:
        raise ValueError('Changed fixed choices or scientific scope')
    if record.get('remaining_obligations')!=REMAINING:
        raise ValueError('An outstanding obligation was removed')
    names={'exact_exterior_defects','signed_edit_ledger','integrable_renormalized_tail',
           'physical_signed_scaling','labelled_heat_and_Fubini_diagnostics'}
    if set(record.get('reviews',{}))!=names or set(record.get('gates',{}))!=names|{'frozen_provenance_preserved','exact_root_dependency_acceptance'}:
        raise ValueError('A review or gate was removed')
    if not all(record['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in record['reviews'].values()):
        raise ValueError('A global moment/stress check failed')
    if not record.get('inherited_provenance') or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited bytes changed')


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
        if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-55):
            raise ValueError('Changed numerical diagnostic at '+path)
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
    print('Exact completed roots and compact inputs remain hypotheses; full pulses, supports and residual remain open.')


if __name__=='__main__':
    main()
