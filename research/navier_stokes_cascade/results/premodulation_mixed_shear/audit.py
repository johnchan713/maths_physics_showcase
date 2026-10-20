#!/usr/bin/env python3
"""Audit actual pre-modulation mixed jets and the retained variance-shape gap."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from mixed_bounds import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='actual-premodulation-mixed-inputs-and-retained-variance-cutoff-budget-gap'
DEPENDENCIES=('fifth_moment_volterra','higher_jet_budget','prescribed_continuation_foundation',
    'annular_exterior_foundation','cutoff_repair_derivative_budget','loop_modulation_bounds')
FLAGS=dict(actual_PREMODULATION_mixed_YE_le_Ahi2_YU_le_Ahi_derived_relative_to_retained_regional_bounds=True,
    native_a_bs_ts_vs_C4_polynomial_bounds_derived_without_pure_fifth_field_input=True,
    earlier_relative_E_bumps_width_point08_C4_roots_and_all_regional_norm_sums_retained=True,
    variance_cutoff_shape_budget_not_supplied_by_manuscript_plateau_clause_identified=True,
    old_loop_H_estimates_require_actual_cutoff_budget_or_separate_composition_estimate=True,
    compressed_family_proves_no_uniform_shape_upper_bound_from_window_even_if_monotone=True,
    shape_dependent_fourth_target_and_unchanged_G_gate_derived=True,
    roots_cutoffs_widths_frequency_supports_gauges_pressure_memory_and_historical_evidence_preserved=True,
    actual_variance_cutoff_shape_identified_with_explicit_radial_sigma_or_Z_budget_verified=False,
    actual_fourth_variance_target_or_old_loop_H_budget_verified_without_shape_premise=False,
    actual_pure_fifth_native_fields_entering_fifth_moments_or_compact_gate_closed=False,
    FINAL_modulated_mixed_jets_or_global_leading_C4_first_C3_second_H23_closed=False,
    full_physical_pulse_covariance_axial_pressure_means_viscous_cutoff_interaction_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=['identify_or_bound_unchanged_variance_cutoff_shape_at_needed_orders_or_prove_actual_composition_budget_without_it',
    'quantify_actual_pure_fifth_native_E_U_and_entering_fifth_moment_data_axis_fifth_bound_already_retained',
    'propagate_pressure_q5_signed_root_and_complete_moving_phase_constants_into_same_compact_gate',
    'complete_FINAL_global_leading_C4_first_C3_second_H23_actual_cutoff_and_retained_scale_inputs',
    'prove_all_physical_pulse_jets_supports_complete_covariance_axial_pressure_means_full_viscous_cutoff_interaction_residual_and_independent_foundations']
def protocol():
    return dict(parent_commit='04a9058da056821a2f4833c10550ec958a9424c1',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_correction='printedp160 C9 specifies smooth bounded variance cutoff and plateau windows,no formula or derivative cap;not the explicit radial sigma;old parameter-only loop jets require shape budget or independent actual composition bound',
        native='same PREMODULATION rectangle Xa exp(tc/2)..rhop exp5,eta[-1,1];Bn=Ahi=A256,A=C4096;originalC4 field norm<=Bn and reference/annularC4 premises retained',
        continuation='Ey=E(1-a)/2,Uy=-kappa beta X nsr/2;referencep1r C4<Bold/100,nsr<=1e8Kold8;no kappa_y for one defining radial slope;referencephysicala=0 afterflattening not substituted foractualprescription',
        annular='earlierrelative E0(1+h),h3coeffs and2additiveU coeffs;d4=5.72e-14,beta=sigma_prime((t-c)/.08+.5),no unit-mass prefactor;beta0<=64,beta_y<=250000;h0=3d4beta0,h1=3d4beta_y;Ey/P<1,Uy<1,aC4<1',
        cover='five regional log-slope C4 bounds Bold,1,.1,1,1000 and Uy bounds Bold,0,1,1,8;SUM not max;log-slope<=Bold+1003<C,Uy<=Bold+10<C,a<=2Bold+2007<C;YE<=Bn2,YU<=Bn',
        algebra='nativea C4<=A,inversea<=A10,bs<=A1286,ts<=A1296,vs<=A2583;positivefloors and complete bs square retained;not FINAL modulated derivatives',
        shape='actual K(t)=zetaL(2+deltaL/8+(deltaL/8)t),Z=max(1,raw shape derivatives1..4);family1-sigma(n(t-.5)+.5) gives H1>=8n,H4>=7424n4 via sigma_third(.5)=-3712 and MVT;not replacement or actual profile failure',
        target='s=zetaL(vs)*sqrt((2+deltaL/2-vs)/a),support deltaL/4<=radicand numerator<=3;C4<=1e12 Z deltaL^-8 A20673;no sqrt of negative value off support',
        gate='logZ+12log10+20673logA+8A20<=A32 sufficient for sC4<=G;Z<=1e100 is auxiliary sufficient example ONLY;1e120A30000deltaL^-8<G for allA>=2^16;unchangedN',
        scope='actual pure fifth fields/entering data variance shape pressure fourth scalar moving phase compact gate FINAL global jets pulse covariance residual and independentfoundations remainopen')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('mixed_shear_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','mixed_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','higher_jet_budget/higher_bounds.py',
        'axis_core_attachment/bounds.py','axis_matching_audit/README.md','compact_jet_envelope/README.md',
        'outer_pressure_pilot/schedule.py','stress_realization_audit/README.md','cutoff_repair_derivative_budget/derivatives.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Mixed shear record changed')
    if not all(expected['gates'].values()):raise ValueError('Mixed shear check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual variance shape, pure fifth fields and FINAL compact gate remain unverified.')
if __name__=='__main__':main()
