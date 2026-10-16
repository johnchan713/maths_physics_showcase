#!/usr/bin/env python3
"""Replay the fixed-width continuation/taper margin checkpoint."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import mpmath as mp
from continuation import budget_review, ode_review, source_review, taper_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='reduced-continuation-relaxed-margins-audited-for-reference-datum'
FLAGS=dict(reference_positive_sources_and_two_sector_growth_derived=True,
           actual_reduced_continuation_and_two_taper_relaxed_margins_derived=True,
           final_constant_shear_segment_closed_by_source_barrier=True,
           complete_forward_pressure_terms_retained=True,
           width_small_field_error_extended_after_final_taper=False,
           whole_compact_shape_join_outer_foundation_independently_verified=False,
           completed_global_backward_pressure_normalization_verified=False,
           unconditional_uniform_actual_pulse_estimate_verified=False,
           unconditional_complete_residual_improvement_verified=False,
           blowup_verified=False)
REMAINING=['audit_shape_transition_join_and_outer_regional_inputs',
           'identify_completed_global_backward_pressure_and_all_exact_moment_identities',
           'verify_coefficient_and_reserved_mean_supports',
           'verify_uniform_actual_pulse_values_derivatives_and_supports_at_all_physical_scales',
           'retain_axial_velocity_and_pressure_in_complete_covariance_mean_corrections',
           'prove_complete_residual_gain_with_viscosity_cutoffs_and_all_interactions']


def protocol():
    return dict(parent_commit='ed1bafb7a24150b8e9353ea25d7c8487e29bba10',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['4.9','4.16','B.22-B.29'],
                fixed_reference_datum_and_parameters_unchanged=True,
                width_and_kappa0='C^-20',Xcut=100,Xi=110,
                continuation_interval='[Xa*exp(C^-20),110]',
                actual_a_lower='C^-22',actual_v_upper='4/5',actual_Pc_lower='11/5',
                reference_source_lower='2.4+.94*Lambda*L*chi',
                reference_growth_range=['2.3','C^3/2'],
                reference_axial_source_error='100*uerr+5*C^-2*exp(B)<j/100',
                comparison_constant='C^4*exp(100*B)',
                comparison_only_through='100*exp(2*C^-20)',
                final_barrier_p1='5/2',final_source_lower='1',
                pressure_and_pressure_eta_kept=True,
                interval_scales=['P16','Md4','Md64'],interval_digits=[100,160],
                log_sum_negligible_tail_interval=['0','1e-200'],
                wide_log_sum_enclosure='max(lower)..max(upper)+log(n)',
                diagnostic_power_exponents=['.5','.1'],diagnostic_digits=75,
                giant_actual_profile_numerically_materialized=False,
                numerical_checks_are_not_continuum_supremum_proofs=True)


def inherited_hashes():
    parent=HERE.parent/'compact_pressure_bridge'
    spec=importlib.util.spec_from_file_location('reduced_frozen_bridge_audit',parent/'audit.py')
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(parent))
    try:
        spec.loader.exec_module(module)
        checks=module.inherited_hashes()
    finally:
        sys.path.pop(0)
    record=json.loads((parent/'evidence.json').read_text())
    for path,expected in record['source_hashes'].items():
        checks['compact_pressure_bridge:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==expected
    return checks


def source_hashes():
    files=[HERE/n for n in ('README.md','continuation.py','audit.py','test_audit.py','protocol.json','.gitignore','requirements.txt')]
    files += [HERE.parent/n/'evidence.json' for n in ('compact_pressure_bridge','core_growth_activation','axis_core_attachment','compact_jet_envelope')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():
        raise ValueError('Frozen protocol changed')
    reviews=dict(exact_taper_identities=taper_review(),positive_source_ODE_comparisons=ode_review(),
                 original_pressure_source_and_final_segment=source_review(),outward_frozen_budgets=budget_review())
    inherited=inherited_hashes()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews.items()}
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
                inherited_provenance=inherited,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def validate(record):
    if record.get('status')!=STATUS or record.get('protocol')!=protocol() or record.get('flags')!=FLAGS:
        raise ValueError('Scope or fixed parameters changed')
    if record.get('remaining_obligations')!=REMAINING:
        raise ValueError('An outstanding obligation was removed')
    names={'exact_taper_identities','positive_source_ODE_comparisons',
           'original_pressure_source_and_final_segment','outward_frozen_budgets'}
    if set(record.get('reviews',{}))!=names or set(record.get('gates',{}))!=names|{'frozen_provenance_preserved'}:
        raise ValueError('A review or gate is missing')
    if not all(record['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in record['reviews'].values()):
        raise ValueError('A continuation gate failed')
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
    result=generate()
    validate(result)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text())
        validate(frozen)
        compare_records(result,frozen)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(STATUS)
    print(f"{len(result['gates'])} gates pass; {len(result['inherited_provenance'])} inherited hashes agree.")
    print('Continuation margins derived for the fixed reference datum; global pressure, pulses and full residual remain open.')


if __name__=='__main__':
    main()
