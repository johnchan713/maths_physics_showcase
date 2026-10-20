#!/usr/bin/env python3
"""Reproduce nested second-coefficient axis jets and the conditional inner join."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from axis_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='nested-complex-core-second-axis-third-jets-and-conditional-inner-collar-join'
DEPENDENCIES=('second_coefficient_source','axis_first_coefficient','core_all_order_foundation',
    'background_jet_foundation','coefficient_support_audit')
FLAGS=dict(strictly_nested_complex_disks_supply_first_raw_jets_for_full_second_source=True,
    leading_third_angular_jet_and_raw_first_second_angular_factorial_retained=True,
    second_core_Picard_bound_retains_true_nested_source_dependent_exponential=True,
    literal_order_two_full_pressure_matrices_and_source_parity_descend_to_X=True,
    conditional_second_axis_safe_X3_Ceta1_estimate_derived=True,
    unchanged_cutoff_equals_one_on_core_and_disjoint_repair_identification_explicit=True,
    conditional_whole_inner_pre_cut_second_third_jet_with_separate_core_collar_bounds=True,
    complete_future_pressure_viscosity_moments_and_all_construction_choices_preserved=True,
    actual_unchanged_cutoff_jets_and_collar_T_evaluated=False,
    finalized_global_second_H23_and_global_repair_angular_inputs_verified=False,
    radial_analyticity_across_flat_activation_or_cutoff_assumed=False,
    all_physical_pulse_derivatives_supports_and_retained_scale_acceptance_closed=False,
    complete_covariance_axial_pressure_means_and_residual_improvement_proved=False,
    independent_unconditional_core_Picard_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['quantify_unchanged_cutoff_jets_and_actual_collar_T_without_replacing_the_cutoff',
    'bound_global_first_U_third_angular_and_leading_U_fourth_angular_inputs_with_complete_second_repair_moments_through_Xv',
    'combine_inner_axis_collar_cutoff_and_global_exact_repair_for_finalized_H23_and_retained_scale_acceptance',
    'prove_all_physical_pulse_jets_supports_full_wave_covariance_axial_pressure_means_and_complete_viscous_cutoff_interaction_residual',
    'independently_check_retained_infinite_core_and_sparse_Picard_regional_hypotheses']
def protocol():
    return dict(parent_commit='f2acf72d1abd5d68653ad25cddbe58507aa90040',
        nesting='first seed |Lambda X|<5;second source/solution |Lambda X|<9/2;first jet disks1/(4Lambda),4.75<5;second core jet disks1/(4Lambda),4.25<4.5',
        inputs='same S>=2,rho<=1,Lambda>0,h<=.01,delta=rho/(800S)<=1/1600;leading B3=2(256S+9+S/Lambda)(1+Lambda)^3;H0=B3 sum_(k=0)^3 k!delta^-k',
        first='same first M1core on3delta/2;H1=6(1+4Lambda)^3(1+2/delta+8/delta^2)*M1core ondelta/2;raw factorial2 retained',
        source='parent full order2 operator/source positive map on complex disk,x=9/(2Lambda),b=1+x;Tcore=max(1,H0,H1),Q<=qTcore^2;current4h lower2h;all pressure/ordered/viscous terms retained',
        inverse='same regular complex ray Gg=xi int_0^1 t^D0 g(txi)dt;D0=(0,0,2,0,3,1);variable diagonal kernel/right derivative sparse annihilation;conditional normal convergence and uniqueness',
        exponential='delta2=delta/4,outer2delta2=delta/2,z=2bqTcore^2,B=6z^2/delta2,M2core=(z+B)exp(B) on3delta2/2;Tcore includes actual first M1core;not set to1',
        parity='literal A0,2 and A1 reflection withP0=(1,1,1,1,-1,-1);full f2 reflection;only first four components even and descend toX',
        core_identity='actual leading seed and finalized first core agree onreal X<=Xa;unchanged kappa=1 throughXa exp(tc/8)>Xa andI3 repair later;forward axis pressure datum never reset',
        join='r<=3:J2_r<=(1+2/delta2)r!(4Lambda)^r M2core+Tcollar^p(A_r M2collar+B_r),p=(0,0,2,4);a=min(1,Xa),a^2<=Xa;x>=4.1/Lambda enclosesXa exp(tc/4) fortc<=.01,b^2>=x;whole pre-cut inner interval and eta stripdelta2/2;separate source/exponential bounds',
        cutoff='same eta-independent second cutoff binomial product with actual k0..k3;k3J0 persists;no global angular repair or pressure moment truncation',
        diagnostics='rational geometric/jet controls,literal symbolic parity/kernels,nine complex manufactured full-source samples;cancel removable symbolic averages beforeX0;axis P2=44/119 matches original hierarchy;not actual PDE or interval proof',
        scope='conditional axis and pre-cut inner second jets only;actual cutoff/collarT,global repair/finalH23,scale,pulses,covariance means,residual and independent foundations open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('axis_second_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','axis_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','axis_first_coefficient/axis_bounds.py',
        'second_coefficient_source/source_bounds.py','core_all_order_foundation/foundation.py',
        'collar_third_coefficient/third_bounds.py','background_summation/hierarchy.py','physical_residual_budget/residual.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Axis second-coefficient record changed')
    if not all(expected['gates'].values()):raise ValueError('Axis second-coefficient check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual cutoff/global repair/H23 and complete pulses remain unverified.')
if __name__=='__main__':main()
