#!/usr/bin/env python3
"""Audit actual PREMOD fifth inputs, same earlier root and fourth pressure."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from native_bounds import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='actual-premod-fifth-fields-memories-and-fourth-pressure-inputs'
DEPENDENCIES=('reference_fifth_transfer','premodulation_mixed_shear','fifth_moment_volterra',
    'prescribed_continuation_foundation','higher_jet_budget','shape_join_pressure',
    'axis_matching_audit','annular_exterior_foundation')
FLAGS=dict(actual_prescribed_continuation_and_shape_C5_derived_from_reference_sources=True,
    same_earlier_eta_dependent_annular_root_C5_le_5p72e14_derived=True,
    whole_PREMOD_pure_fifth_FE_le_A2_FU_le_A_derived_relative_to_retained_core_and_value_lemmas=True,
    genuine_six_entering_fifth_traces_and_full_native_C5_moments_bounded_by_A8=True,
    full_native_fourth_pressure_rows_bounded_by_distinct_A1546_envelope=True,
    axis_radial_pressure_U_squares_normalizers_and_all_regional_suprema_retained=True,
    no_historical_record_root_width_frequency_support_or_pressure_gauge_changed=True,
    actual_variance_shape_budget_or_old_parameter_only_loop_H_budget_accepted=False,
    fourth_scalar_q5_signed_root_moving_phase_or_compact_gate_closed=False,
    FINAL_modulated_repaired_C4_C5_global_first_C3_second_H23_or_all_scale_pulses_closed=False,
    complete_covariance_axial_pressure_means_full_viscous_cutoff_interaction_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['propagate_actual_native_pressure_A1546_and_retained_mixed_shears_through_fourth_signed_root_scalar_q5_and_complete_moving_phase_maps',
    'bound_unchanged_variance_cutoff_shape_or_prove_independent_actual_composition_budget_before_parameter_only_loop_H_or_compact_acceptance',
    'close_unchanged_compact_gate_then_FINAL_global_leading_C4_first_C3_second_H23_actual_cutoffs_and_retained_scale_inputs',
    'prove_all_scale_physical_pulse_jets_common_supports_complete_covariance_axial_pressure_means_full_viscous_cutoff_interaction_residual_and_independent_foundations']
def protocol():
    return dict(parent_commit='0963f1379a76c213f16db60a5acb52c3cd0cb658',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        domain='actual PREMOD native input Xminus=Xa exp(tc/2)..Xrep=rhop exp5,eta[-1,1];fields extend to axis only for forward moments;NOT FINAL modulation or compact repair',
        parameters='Kold,Bold,Xa4/Lambda,Xi110,deltaC^-20,Tsh2000B,XR110(CP)^10,Xsep110expTsh,A=C4096>=2^16,Bn=Ahi=A256;all choices unchanged',
        continuation='same actual prescribed slopes including kappa beta theta and final .8 shear;logCE C5<=1e11LambdaK25<B/100;gC5<=j+K/Lambda+1e8K3delta<3e-18;no inverse kappa or reference zero shear shortcut',
        moments='core CE/sqrtX C5<expB,shape Eb=exp204B/C;all five genuine axis-forward moments and same Pi_axis retained;Cp core integrable then finite log length;Cp C5<=C^-2exp412B',
        annular='C5 inversef5,inversef2=25,energy8125;target<2.86e-17;eta-dependent quadratic C5<=100,inverse1000,same rootC5<=5.72e-14;relative E bumps width.08 and value64 unchanged',
        cover='SUM continuation_shape,ideal_restoration,annular,outerprefix bounds;E5<=1+21P+22P+1000A<A2,U5<=36<A;FE<=A2,FU<=A;all six actual fifth entrance traces<=A8',
        pressure='retained C4 Bn and YE<=Bn2,YU<=Bn;complete Volterra B1,B2<=A517;R0=10000A1544;native p1,p2 C4<=A1546,NOT assigned old Ahi=A256;axis Cp and U-square memory retained',
        scope='conditional analytical derivations with finite manufactured rational and two-precision interval checks;variance cutoff budget,fourth scalar/q5/moving phase,compact acceptance,FINAL global jets,pulses,covariance,residual and independent foundations remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('actual_fifth_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','native_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','reference_fifth_transfer/reference_bounds.py',
        'higher_jet_budget/higher_bounds.py','axis_core_attachment/bounds.py',
        'prescribed_continuation_foundation/continuation_checks.py','core_all_order_foundation/coordinate_checks.py',
        'compact_jet_envelope/README.md','fifth_moment_volterra/volterra.py')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
@lru_cache(None)
def reviews():return review()
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g=dict(r['checks']);g['historical_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Actual fifth record changed')
    if not all(expected['gates'].values()):raise ValueError('Actual fifth check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual variance shape, fourth moving maps and FINAL compact gate remain open.')
if __name__=='__main__':main()
