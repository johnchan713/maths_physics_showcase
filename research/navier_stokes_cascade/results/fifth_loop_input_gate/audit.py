#!/usr/bin/env python3
"""Tame fifth moving loop: explicit actual-native and raw-cutoff gate."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from moving5 import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='tame-fifth-moving-loop-original-native-and-cutoff-input-gate'
PARENT='0657b89fddc9b290f204f3d0fd719238949ad61f'
DEPENDENCIES=('final_fifth_moment_interface','joint_fourth_loop_gamma',
    'fourth_loop_shape_gate','native_fifth_mixed_cover','native_mixed_fourth_boundary',
    'actual_cutoff_pullback','modulation_repair_all_order','repair_state_bounds')
MANIFEST_SHA='076c320e38144a3d00115a13dc84451007039990e187633f872e6c7f9137eab0'
FLAGS=dict(tame_fifth_signed_root_and_full_moving_inverse_affine_in_top_inputs_derived=True,
    q6_required_by_fifth_removable_R_and_signed_origin_retained=True,
    original_native_and_raw_cutoff_sufficient_B5_inequality_derived=True,
    native_total_order_six_and_full_pressure_sixth_principal_rows_identified=True,
    full_incoming_axis_Cp_and_nonzero_axial_U_memory_retained=True,
    same_C_K_N_widths_roots_cutoffs_supports_gauges_history_unchanged=True,
    actual_R5_Z5_T5_Lambda5_or_B5_enclosed=False,
    actual_joint_fourth_target_J_or_raw_cutoff_budget_verified=False,
    fifth_y_or_full_joint_three_variable_fifth_loop_bound_derived=False,
    unconditional_compact_or_global_leading_C4_accepted=False,
    uniform_actual_pulses_complete_covariance_axial_pressure_means_residual_or_blowup_verified=False,
    independent_foundation_review_or_formal_certificate=False)
REMAINING=[
    'enclose_actual_native_R5_from_E_U_1_5_and_full_M_S_Pi_0_6_or_direct_scalar_propagation;and_actual_Z5_or_direct_T5_Lambda5;B5_not_verified',
    'independently_verify_joint_fourth_target_J_or_retained_raw_cutoff_gate;B5_does_not_replace_J;then_join_complete_compact_maps_to_global_regions_and_physical_scales',
    'prove_coefficient_reserved_mean_supports_uniform_actual_pulses_complete_covariance_axial_pressure_means_and_full_viscous_cutoff_interaction_residual']

def protocol():
    return dict(parent_commit=PARENT,
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        premises='retained native/core/regional/floor/value/root premises AND unverified parent J:target jointJ4<=G;lower scalar_nativeC4<G,muC4<exp600G,ThetaJ4<exp80000G;no final theorem premise',
        inputs='R5=1+sup|a5|+sup|bs5|+sup|E5|+sup|p2_5|;T5=R5+sup|s5|;raw Z5=max(1,sup|Kshape^(j)|,1..5),not factorial;no actual enclosure',
        signed_root='qJ6<exp4G,RJ5<exp5G,R2C5<exp10G,sqrtC5<exp71G,FpartialC5<exp100G,Fmu_inverse<exp10G;mu5=(s5-Fp*p5-lower)/Fmu<exp3200G*T5;no p,mu,Vmu division',
        moving='degree-five eta_theta or eta_externalphase jets at FIXED y,uniform in y;top jets linear;fixed_t5<exp20000G*T5,fixed_sources_w5<exp100000G*T5,Phi5<exp100010G*T5,positiveTheta5<exp600010G*T5,movedsources5<exp700000G*T5;no new fifth_y rows',
        primitive='exact weighted zero means with E_bs and variance;fixed_external_phase zero_mean_integral cost2_each,sum4;Lambda5<exp1000000G*T5;not physical covariance mean closure',
        cutoff='active radicand>=deltaL/4;scalar vpartialC5<=100G8,raw_cutoff_sqrt_product<1e10 Z5 deltaL^-10;target partialC5<1e24G48 Z5 deltaL^-10<expG Z5 deltaL^-10;tame s5<exp1000G Z5 deltaL^-10 R5',
        gate='Lambda5<exp(1100000G+10A20)*Z5*R5;B5:logZ5+logR5+10A20+1100000G<=G2;J_AND_B5=>Lambda5<expG2<H=>same_F5=>compactordinaryC4<H925<H1024;not K replacement,global acceptance or J verification',
        native_boundary='a5=-2E_(1,5)/E+lower,bs5=2U_(1,5)/E+lower;p2_5=6/(LE)*[((1-eta2)U-Deta)M6-(1-eta2)S6-X(1-eta2)Pi6]+lower;Pi6 full_axis_plus_Cp,no axis reset;nonzero_U*M6 retained,p1 not signed root argument',
        scope='conditional tame original-input gate ONLY;actual R5_Z5_T5_Lambda5_B5_and_J remain unverified;history unchanged;finite exact/outward diagnostics not continuum proof or actual profile failure;no global pulse covariance residual blowup')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m

@lru_cache(None)
def history():
    p=HERE.parent/'final_fifth_moment_interface';prior=sys.path[:];sys.path.insert(0,str(p))
    try:parent=module('fifth_input_gate_parent_audit',p/'audit.py')
    finally:sys.path[:]=prior
    names=sorted(parent.HISTORICAL_RECORDS+('final_fifth_moment_interface',))
    digest=hashlib.sha256(json.dumps(names,separators=(',',':')).encode()).hexdigest()
    if len(names)!=135 or len(set(names))!=135 or digest!=MANIFEST_SHA:raise ValueError('Fixed 135-record parent manifest changed')
    p=HERE.parent/'harmonic_closure_foundation'
    prior=sys.path[:];sys.path.insert(0,str(p))
    try:m=module('fifth_input_gate_historical_binding_checker',p/'audit.py')
    finally:sys.path[:]=prior
    h=m.verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
    bindings=h['bindings']
    if not bindings or not all(bindings.values()):raise ValueError('Historical source binding failed')
    return dict(all_bindings_verified=True,binding_count=len(bindings),
        binding_result_sha256=hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest(),
        fixed_manifest_sha256=digest,
        frozen_records={str((HERE.parent/n/'evidence.json').relative_to(PROJECT)):sha(HERE.parent/n/'evidence.json') for n in names})

def sources():
    files=[HERE/n for n in ('README.md','moving5.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py',
        'final_fifth_moment_interface/audit.py','final_fifth_moment_interface/fifth.py',
        'joint_fourth_loop_gamma/joint.py','fourth_loop_shape_gate/loop_fourth.py')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in files}

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(review());h=history();g=dict(r['checks']);g['historical_source_bytes']=h['all_bindings_verified']
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Fifth loop original-input record changed')
    if not all(expected['gates'].values()):raise ValueError('Fifth loop original-input check failed')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual native/cutoff B5 and J remain unverified.')
if __name__=='__main__':main()
