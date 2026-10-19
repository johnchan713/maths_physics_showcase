#!/usr/bin/env python3
"""Reproduce sharp smooth-shape feasibility, not acceptance of the actual cutoff."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from shape_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='sharp-smooth-cutoff-capacity-and-eliminated-affine-budget-test'
DEPENDENCIES=('cutoff_cone_scale_gate','cutoff_repair_derivative_budget','first_seminorm_bound',
    'global_swirl_envelope','background_gap_certificate','leading_cone_margin')
FLAGS=dict(sharp_joint_derivative_enclosure_capacity_derived_for_the_entire_smooth_cutoff_class=True,
    necessity_does_not_assume_monotonicity_and_strictness_uses_flat_endpoint=True,
    auxiliary_smooth_existence_established_by_interior_shrink_mollify_and_normalize=True,
    sharp_weighted_viscous_row_infimum_is_strict_and_unattained=True,
    exact_rational_quadratic_eliminates_shape_budgets_for_the_sufficient_affine_certificate=True,
    upper_enclosures_are_not_confused_with_exact_attained_norms=True,
    complete_corrected_width_pressure_memory_and_viscosity_rows_retained=True,
    original_cutoff_coefficients_profiles_roots_frequency_tilt_scale_gauges_and_mean_memories_unchanged=True,
    historical_frozen_evidence_and_original_comparators_preserved=True,
    auxiliary_feasible_shape_or_budget_is_the_unchanged_actual_cutoff=False,
    impossible_sufficient_bound_proves_actual_cone_or_PDE_failure=False,
    actual_cutoff_shape_existing_tilt_or_retained_scale_acceptance_verified=False,
    weighted_stress_row_infimum_is_necessarily_infimum_of_full_max_of_rows=False,
    actual_huge_profile_or_physical_pulses_numerically_resolved=False,
    retained_Navier_Stokes_core_and_regional_premises_independently_proved_unconditional=False,
    unconditional_pulse_complete_residual_improvement_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['identify_valid_upper_enclosures_of_the_unchanged_actual_cutoff_not_an_auxiliary_shape_then_check_the_sharp_joint_capacity',
    'verify_same_existing_tilt_margin_and_retained_scale_with_full_rows_correct_tc_and_all_other_acceptance_gates',
    'prove_uniform_actual_physical_pulse_size_derivative_support_and_inverse_gap_on_every_scale',
    'retain_complete_wave_covariance_axial_velocity_pressure_and_means_in_full_viscous_cutoff_interaction_residual',
    'independently_verify_retained_exact_core_and_regional_foundations_before_unconditional_claims']
def protocol():
    return dict(parent_commit='355dce1e76bf0e6d97ec7da84aa06e0fddff5522',
        class_specification='smooth K=1 on s<=0,K=0 on s>=1,0<=K<=1; M,N are upper enclosures of sup|Kprime|,sup|Ksecond|, not necessarily exact norms',
        retained_geometry='s=log(X/Xkeep)/Delta,Delta=tc/16=C^-20/(160 sqrt(log C)); allowed window and actual cutoff unchanged',
        capacity='A=integral_0^1 min(M,Ns,N(1-s));A=M-M^2/N when 2M<=N,else N/4; existence of SOME smooth shape iff A>1',
        smooth_witness='eps=(A-1)/(16A);shrink height and support by 1-2eps,mass m=(1-2eps)^2A>1; nonnegative smooth convolution supported within eps/2; normalize mass; upper jets (1-2eps)M/m,N/m',
        weighted_cost='Rw>A0+A1+2A2+2sqrt(A2(A1+A2)); sharp unattained infimum for the stress row only, not necessarily the full max',
        eliminated_gate='constant rows<=B;cap=min(1,(B-D0)/D1-1)>0;beta=B-A0-A1-2A2>0;P=A1+A2,Q=A2;t0=min(cap,beta/(2P));P*t0^2-beta*t0+Q<0',
        rational_budget_witness='M=1+t0;Nlo=M^2/t0;Nhi=(B-A0-A1M)/A2;N=(Nlo+Nhi)/2; only auxiliary shape enclosure existence',
        negative_controls='old relaxed corner can pass while no smooth shape satisfies the certificate; exact capacity/frontier equality rejected; loose upper enclosures need not obey N>=M^2',
        diagnostics='25 exact capacity controls,4 strict frontier/padding controls,15 weighted infimum approximations,symbolic identities and complete inherited row composition; not actual-profile samples',
        scope='new continuum cutoff lemma; its Navier-Stokes application conditional on retained full row estimates; actual cutoff,tilt,scale,pulses,complete residual and blowup not verified')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('smooth_shape_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','shape_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','cutoff_cone_scale_gate/gate_bounds.py',
        'cutoff_repair_derivative_budget/derivatives.py','first_seminorm_bound/seminorms.py','annular_first_correction/README.md','core_growth_activation/README.md')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
@lru_cache(None)
def reviews():return review()
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g=dict(r['checks'])
    g['historical_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Smooth cutoff feasibility record changed')
    if not all(expected['gates'].values()):raise ValueError('Smooth cutoff feasibility check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; auxiliary shape feasibility does not verify the unchanged cutoff or scale.')
if __name__=='__main__':main()
