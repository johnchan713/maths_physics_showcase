#!/usr/bin/env python3
"""Replay the actual shape/join-entry and exact leading pressure ledger."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import mpmath as mp
from shape import budget_review, jet_review, pressure_review, scalar_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='shape-join-entry-and-exact-leading-pressure-ledger-audited'
FLAGS=dict(actual_shape_and_prejoin_relaxed_margins_derived=True,
           actual_C2_incoming_join_entry_derived=True,
           same_exact_C2_join_root_accepted=True,
           leading_backward_pressure_identity_for_exact_repair_roots_derived=True,
           exact_outer_Cp_correction_maps_are_explicit_dependencies=True,
           all_global_moment_and_signed_stress_identities_verified=False,
           whole_leading_profile_independently_verified=False,
           coefficient_and_reserved_mean_supports_verified=False,
           unconditional_uniform_actual_pulse_estimate_verified=False,
           unconditional_complete_residual_improvement_verified=False,
           blowup_verified=False)
REMAINING=['join_outer_regional_cone_inputs_to_inner_shape_join_certificate',
           'audit_renormalized_I_and_global_M_J_S_and_signed_stress_identities',
           'verify_coefficient_and_reserved_mean_supports',
           'verify_uniform_actual_pulse_values_derivatives_and_supports_at_all_physical_scales',
           'retain_axial_velocity_and_pressure_in_complete_covariance_mean_corrections',
           'prove_complete_residual_gain_with_viscosity_cutoffs_and_all_interactions']


def protocol():
    return dict(parent_commit='a8371b14e525bbb03656909f90120221ae63009f',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['4.15-4.16','A.11','B.32-B.39','C.17'],
                same_fixed_parameters_datum_and_exact_roots=True,
                shape_length='2000*B',shape_source_lower='1',shape_Pc_lower='5/2',
                shape_a_range=['.7','.9'],shape_vs_upper='.9',
                entry_C2_upper='.255*epsilon',target_C2_upper='.286*epsilon',
                exact_coefficient_C2_upper='5.72e-14',
                inverse_f_squared_C2_norm=20,raw_second_derivative_factor=2,
                join_radius_floor=10000,join_Pc_lower='3.32',join_vs_upper='1',
                backward_pressure_scope='completed leading swirl using defined exact repair roots, not full corrected pressure',
                Cp_edit_totals=['angular:0','axial:0','heat+compensation:0','axis+join:0','modulation+repair:0'],
                pressure_equality_is_functional_in_eta=True,
                compact_modulation_input_hypotheses_remain_dependencies=True,
                raw_I_infinite_tail='diverges as X^(1-h); use renormalized I',
                interval_scales=['P16','Md4','Md64'],interval_digits=[100,160],
                manufactured_fixtures_are_actual_profile=False,
                diagnostics_are_continuum_supremum_certificates=False)


def inherited_hashes():
    parent=HERE.parent/'reduced_continuation_margins'
    spec=importlib.util.spec_from_file_location('shape_frozen_continuation_audit',parent/'audit.py')
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(parent))
    try:
        spec.loader.exec_module(module)
        checks=module.inherited_hashes()
    finally:
        sys.path.pop(0)
    record=json.loads((parent/'evidence.json').read_text())
    for path,expected in record['source_hashes'].items():
        checks['reduced_continuation_margins:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==expected
    return checks


def source_hashes():
    files=[HERE/n for n in ('README.md','shape.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/'evidence.json' for n in ('reduced_continuation_margins','axis_matching_audit','angular_c2_transfer','pulse_moment_audit','heat_exterior_audit','repair_state_bounds')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():
        raise ValueError('Frozen protocol changed')
    reviews=dict(shape_source_and_C2_root_scalars=scalar_review(),
                 complete_C2_normalization_and_energy_cancellation=jet_review(),
                 exact_pressure_ledger_and_renormalization_controls=pressure_review(),
                 outward_actual_shape_join_and_Cp_root_budgets=budget_review())
    inherited=inherited_hashes()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews.items()}
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
                inherited_provenance=inherited,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def validate(record):
    if record.get('status')!=STATUS or record.get('protocol')!=protocol() or record.get('flags')!=FLAGS:
        raise ValueError('Changed fixed choices or scientific scope')
    if record.get('remaining_obligations')!=REMAINING:
        raise ValueError('An outstanding obligation was removed')
    names={'shape_source_and_C2_root_scalars','complete_C2_normalization_and_energy_cancellation',
           'exact_pressure_ledger_and_renormalization_controls','outward_actual_shape_join_and_Cp_root_budgets'}
    if set(record.get('reviews',{}))!=names or set(record.get('gates',{}))!=names|{'frozen_provenance_preserved'}:
        raise ValueError('A review or gate was removed')
    if not all(record['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in record['reviews'].values()):
        raise ValueError('A shape/join/pressure gate failed')
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
        if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-70):
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
    print('Leading pressure uses exact roots; outer/global stress identities, physical pulses and full residual remain open.')


if __name__=='__main__':
    main()
