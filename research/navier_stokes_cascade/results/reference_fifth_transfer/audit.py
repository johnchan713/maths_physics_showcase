#!/usr/bin/env python3
"""Audit sixth reference pressure inputs and fifth integrated source transfer."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from reference_bounds import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='sixth-reference-pressure-and-fifth-source-transfer-with-full-memory'
DEPENDENCIES=('premodulation_mixed_shear','fifth_moment_volterra','higher_jet_budget',
    'prescribed_continuation_foundation','core_all_order_foundation','axis_core_attachment')
FLAGS=dict(axis_positive_mixture_C6_le_1050P2_derived=True,
    same_monotone_reference_taper_preserves_C6_without_seventh_mixed_input=True,
    reference_regular_F_E_U_C6_and_full_pressure_C6_bounded_relative_to_retained_core=True,
    actual_core_entry_Xa_five_fifth_moment_inputs_bounded_without_axis_reset=True,
    reference_Sq_Sn_and_integrated_p1_ns_C5_bounded_in_unchanged_B_budget=True,
    exact_rational_source_ODE_diagnostics_retain_axis_and_radial_pressure=True,
    radial_sigma_not_confused_with_unspecified_variance_cutoff=True,
    roots_supports_widths_frequency_gauges_and_historical_evidence_preserved=True,
    actual_prescribed_continuation_shape_join_and_full_PREMOD_pure_fifth_gate_closed=False,
    actual_variance_shape_or_old_parameter_only_loop_H_budget_verified=False,
    fourth_pressure_q5_moving_phase_and_compact_gate_accepted=False,
    FINAL_global_leading_C4_first_C3_second_H23_or_uniform_physical_pulse_bounds_closed=False,
    complete_covariance_axial_pressure_means_full_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['propagate_reference_C5_sources_through_actual_prescribed_continuation_shape_and_same_annular_root_and_outer_prefix',
    'carry_same_profile_entering_fifth_moments_to_actual_compact_left_endpoint_then_evaluate_full_forward_Volterra_input',
    'bound_unchanged_variance_cutoff_shape_or_independently_bound_actual_composition_without_importing_radial_sigma',
    'quantify_fourth_pressure_q5_signed_root_and_complete_moving_phase_constants_and_unchanged_compact_gate',
    'complete_FINAL_global_jets_uniform_all_scale_pulses_supports_full_covariance_axial_pressure_means_and_viscous_cutoff_interaction_residual']
def protocol():
    return dict(parent_commit='7a282097f781d13dc410e64d606109aff676f6c2',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        domain='reference and native X0=Xa=4/Lambda through Xi=110,eta[-1,1];not actual prescribed profile or final compact modulation',
        norms='factorial angular Ck=sum sup abs derivative_j/j!,K=Kold;0<h<=.01;retained Md64 parameters and radial width C^-20 unchanged',
        axis='exact positive measure mass<5;Pi_axis=-P2 integral f^(2theta),f=1/(1+eta2);coefficient cap binomial(j+3,3);C6<=1050P2',
        reference='B22 uses explicit radial 1-sigma;logRr and Ur are positive averages of native logPhi and Unat;C6 needs no mixed seventh input or derivative of radial cutoff',
        tame='Rr C6<=3K7,inverse<=15K6;regular F C6<=1e8 exp(-Lambda)Lambda6K49<.01;E C6<1;Ur C6<K;complete Pi C6<=1200P2<K',
        sources='Sq C5<=1e4 Lambda K10,Sn C5<=1e3K2;Pi_eta needs Pi C6;axis and cumulative Cp both retained',
        transfer='p1r C5<=1e8 X Lambda K24,nsr C5<=1e5K3,both<Bold/100;variation formula cancels phi_star but not its angular gradient or pressure memory',
        entry='at core Xa only:M C5<=aK,I<=a2,J<=a2K,S<=aK2+a2/2,Cp<=a;all<K;not reset or transferred to later actual entry',
        scope='conditional analytic estimates plus finite exact rational/interval diagnostics;variance shape,actual full fifth fields,compact gate,FINAL jets,pulses,covariance,residual and independent foundations remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('reference_fifth_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','reference_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','higher_jet_budget/higher_bounds.py',
        'axis_core_attachment/bounds.py','prescribed_continuation_foundation/continuation_checks.py',
        'core_all_order_foundation/coordinate_checks.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Reference fifth record changed')
    if not all(expected['gates'].values()):raise ValueError('Reference fifth check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual full fifth fields, variance shape and compact gate remain open.')
if __name__=='__main__':main()
