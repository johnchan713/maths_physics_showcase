#!/usr/bin/env python3
"""Reproduce exact fourth pressure/moving-phase input-boundary findings."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from input_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='exact-fourth-pressure-principal-rows-and-moving-phase-input-boundary'
DEPENDENCIES=('fourth_target_normalization','higher_jet_budget','modulation_repair_all_order',
    'compact_pressure_bridge','prescribed_continuation_foundation')
FLAGS=dict(exact_normalized_fifth_to_fourth_pressure_principal_combinations_derived=True,
    actual_axial_U_and_full_denominator_lower_jets_retained=True,
    complete_incoming_Cp_and_axis_pressure_fifth_jets_retained=True,
    compatible_forward_moment_family_proves_C4_input_to_C4_pressure_nonimplication=True,
    full_fourth_moving_inverse_composition_and_primitive_gate_transfer_derived=True,
    p2_only_variance_input_and_signed_root_zero_limit_distinguished=True,
    roots_width_frequency_supports_gauges_historical_evidence_preserved=True,
    actual_profile_fails_or_counterexample_is_selected_cone_profile=False,
    actual_fifth_moment_combinations_mixed_shear_inputs_or_scalar_fourth_constants_quantified=False,
    actual_a4_b4_compact_gate_global_leading_U_C4_first_U_C3_or_H23_closed=False,
    full_physical_pulse_scale_support_covariance_mean_or_complete_residual_verified=False,
    unconditional_foundations_blowup_or_formal_certificate=False)
REMAINING=['quantify_actual_T2_fifth_moment_axis_pressure_combination_and_mixed_y_eta4_shear_inputs_without_resetting_memory',
    'bound_fourth_signed_root_scalar_partials_including_q_order5_and_all_moving_phase_terms_at_same_parameters',
    'evaluate_LA_LB_and_sufficient_2alphaLA_plus_2betaLB_gate_at_fixed_N_for_actual_global_leading_C4',
    'supply_finalized_first_C3_second_H23_full_pressure_repair_actual_cutoff_and_retained_scale_inputs',
    'prove_all_physical_pulse_derivative_supports_complete_covariance_axial_pressure_means_viscosity_cutoff_interaction_residual_and_independent_foundations']
def protocol():
    return dict(parent_commit='187149d4cdd64da1d82d405dc3cfbfd3c92519f5',
        pressure='original bridge.forward;Pi=Pi_axis+Cp;normalized moment fifth jets m5,i5,j5,s5,c5,k5;R1,R2 evaluate ONLY these fifth coefficients as zero,all lower denominators intact',
        principal='p1_4=R1+5/L*(d*m5-(D*eta*i5+d*j5)/Hx);p2_4=R2+5/(L*E)*((d*U-D*eta)*m5-d*s5-X*d*(c5+k5));factor5=5!/4!',
        counterexample='manufactured patch 1<=X<=2,E=1,U=eps*n^-4*sin(n*(eta-1/2));all exact forward moments and C4 norms uniform;delta(p2)_4/n -> -D*eps/(48*(1-h/2));not selected profile or cone failure',
        root='same signed F(mu,p2)=s;mu4=(s4-Fp*p2_4-RF)/Fmu;Fmu(0)=d0/sqrt2;never Vmu;scalar partials4 require q5;not a new varying p1 argument',
        moving='q4=-(Phi40/24+R4(Phi))/w;g_pullback4=g40/24-g01*Phi40/(24*w)+R4(g)-g01*R4(Phi)/w;all fourth Taylor partitions and complete phase integral retained',
        primitive='period-one zero-mean antiderivative sup norm<=2;actuala4<=2LA,b4<=2LB;LA/LB are full source expressions of equation7,not evaluated',
        gate='parentalpha=3200000*A*Za,beta=1+3200000*A*Zb,R=N*(K-Kknown)-(1+3200000*A*ZH)*H;2alphaLA+2betaLB<=R sufficient only;N and all parameters unchanged',
        checks='35 exact original-forward principal splits,4 exact compatible-field controls and symbolic infinite-n limit,all fourth mixed phase partitions;inverse toy raw405/512 checked by Lagrange inversion differs from previous explicit-lift raw9/512;not continuum proof',
        scope='actual native higher inputs/scalar constants/compact gate globalC4/H23 physicalpulse covariance means residual and independentfoundations remainopen')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('fourth_input_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','input_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','compact_pressure_bridge/bridge.py',
        'fourth_target_normalization/target_bounds.py','higher_jet_budget/higher_bounds.py','compact_jet_envelope/jets.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Fourth input record changed')
    if not all(expected['gates'].values()):raise ValueError('Fourth input check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual higher inputs, compact gate and global H23 remain unverified.')
if __name__=='__main__':main()
