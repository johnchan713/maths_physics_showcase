#!/usr/bin/env python3
"""Complete second inner source and conditional positive-collar third jets."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from source_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='full-second-inner-source-derivative-envelope-and-conditional-positive-collar-J3'
DEPENDENCIES=('axis_first_coefficient','collar_third_coefficient','background_summation',
    'higher_pressure_source_moments','cutoff_higher_jet_interface','collar_operator_envelope')
FLAGS=dict(full_second_radial_pressure_theta_and_axial_sources_with_xi_derivatives_bounded=True,
    leading_extra_third_eta_derivative_and_regular_velocity_cost20_explicit=True,
    lower_lambda1_current_lambda2_and_leading_radial_axial_shift_distinguished=True,
    conditional_positive_collar_second_inner_third_jet_with_true_T_dependent_M2_derived=True,
    whole_inner_first_input_uses_actual_cutoff_jets_and_factorial_to_raw_conversion=True,
    complete_pressure_viscosity_interaction_memory_and_global_repair_obligations_retained=True,
    historical_evidence_construction_gauges_roots_width_frequency_schedule_supports_unchanged=True,
    actual_unchanged_cutoff_jets_and_resulting_T_evaluated=False,
    second_axis_jets_global_repair_or_finalized_H23_verified=False,
    local_source_collar_estimate_truncates_global_pressure_moments=False,
    all_physical_pulse_jets_supports_and_retained_scale_acceptance_closed=False,
    complete_covariance_axial_pressure_means_and_residual_improvement_proved=False,
    independent_unconditional_core_Picard_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['quantify_unchanged_first_cutoff_jets_and_actual_T_from_same_core_collar_inputs',
    'complete_second_axis_cutoff_and_global_repair_with_full_future_pressure_moments_to_actual_H23',
    'quantify_global_angular_repair_inputs_and_verify_retained_value_scale_acceptance',
    'prove_all_physical_pulse_derivative_support_estimates_full_covariance_axial_pressure_means_and_complete_viscous_cutoff_interaction_residual',
    'independently_check_retained_infinite_core_sparse_Picard_and_regional_hypotheses']
def protocol():
    return dict(parent_commit='6bc32141ce0df0ebe51f6d1ea55c599c232166a5',
        source='w1=Omega1/X with both (0,1),(1,0),radial viscosity and -Z_(-1-D)Z_-1 v0;P2=F1^2-w1/2;full theta/axial sources with lower exponents',
        exponents='current matrix lambda2=4h;lower lambda1=2h;radial previous a=-1;theta first -1+h;axial first -1/2+h;current pressure coefficient -1+2h',
        derivatives='leading X<=3 eta<=3;first X<=3 eta<=2;v0 eta2 needs Ubar0 eta3;all needed leading v costs<=20 and first v eta1<=12 by uniform endpoint majorant',
        leading_upgrade='same K3 E on4delta;H0=K3 E sum_(k=0)^3 k! delta^-k on2delta;no changed activation or new radial derivative input',
        first_input='actual eta-independent cutoff binomial jets through3 times whole-inner first Ceta2 bound;multiply by2 for raw eta-second jets;on delta/2;disjoint I3 repair cannot feed backward into forward inner averages',
        pair='positive pair(value,xi gradient),X=(x,2b),xi=(b,1);each full source/row polynomial degree<=2 in H0,H1;T scalar>=max(1,H0,H1),Q2<=qT^2,Q2xi<=qxiT^2',
        solution='retained delta<=1/1600,delta2=delta/4;z2=2bqT^2,B2=6z2^2/delta2;M2=(z2+B2)exp(B2) on3delta2/2;NOT evaluated atT1',
        third='reuse same D0/sparseA1 xi differentiation on xi>=a>0;Ceta1 J2_r<=T^p(A_r M2+B_r),p=(0,0,2,4);A,B calculated atq,qxi,delta2;not a second-axis or finalized H23 bound',
        diagnostics='nine manufactured original hierarchy source/value/gradient comparisons;symbolic U0=eta^3 produces w1=6 at eta0;exponent/viscosity omission controls;not actual PDE fields or interval certification',
        scope='actual cutoff norms,T,second axis/global repair/final H23,physical pulse/scale/covariance/mean/residual and independent foundations remain open;no global pressure truncation or reset')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('second_source_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','source_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','axis_first_coefficient/axis_bounds.py',
        'collar_third_coefficient/third_bounds.py','background_summation/hierarchy.py','physical_residual_budget/residual.py',
        'cutoff_higher_jet_interface/jet_bounds.py','core_all_order_foundation/foundation.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Second-source record changed')
    if not all(expected['gates'].values()):raise ValueError('Second-source check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual cutoff/T, global second H23 and complete pulses remain unverified.')
if __name__=='__main__':main()
