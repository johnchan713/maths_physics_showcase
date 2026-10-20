#!/usr/bin/env python3
"""Reproduce fifth-moment Volterra and fourth pressure/shear input estimates."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from volterra import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='fifth-forward-moment-volterra-kernels-and-fourth-input-majorants'
DEPENDENCIES=('fourth_pressure_input_boundary','fourth_target_normalization','higher_jet_budget',
    'prescribed_continuation_foundation','annular_exterior_foundation')
FLAGS=dict(exact_both_fifth_moment_combinations_reduced_to_pure_field_jets_and_incoming_memory=True,
    energy_and_complete_pressure_combination_gives_quadratic_local_swirl_kernel_bound=True,
    actual_axial_U_square_Cp_and_axis_datum_terms_not_removed=True,
    explicit_positive_fourth_pressure_lower_remainder_and_mixed_shear_majorants_derived=True,
    no_C4_to_C5_or_all_order_smoothness_to_smallness_promotion=True,
    roots_width_frequency_supports_gauges_and_historical_evidence_preserved=True,
    actual_pure_fifth_field_mixed_shear_or_entering_fifth_inputs_quantified=False,
    actual_fourth_scalar_variance_moving_phase_constants_or_compact_gate_accepted=False,
    global_leading_C4_first_C3_second_H23_or_all_physical_pulse_bounds_closed=False,
    full_covariance_axial_pressure_means_viscosity_cutoff_interaction_residual_or_blowup_verified=False,
    formal_certificate=False)
REMAINING=['quantify_actual_FE_FU_YE_YU_and_six_entering_fifth_jets_through_same_native_continuation_shape_and_eta_dependent_annular_root',
    'propagate_explicit_native_majorants_through_fourth_signed_root_scalar_q5_and_complete_moving_phase_maps',
    'evaluate_unchanged_2alphaLA_plus_2betaLB_gate_and_actual_global_leading_C4',
    'finish_first_C3_second_H23_full_pressure_repair_actual_cutoffs_retained_scale_and_all_physical_pulse_support_derivatives',
    'prove_complete_covariance_axial_pressure_means_full_viscous_cutoff_interaction_residual_and_independent_foundations']
def protocol():
    return dict(parent_commit='e85d339f19d696a915a047a923cb0303f4824acd',
        moments='all normalized fifth jets;M,I,J,S,Cp genuine forward integrals from fixed radial a with all incoming fifth traces;Pi_axis retained;Lu,Le,Lue contain every lower convolution',
        T1='Tin1+int d(1-sqrt(2x)E/Hx)u5-int sqrt(2x)(Deta+dU)e5/Hx-d int sqrt(2x)Lue/Hx',
        T2='Tin2+int[d(U(X)-2U(x))-Deta]u5+d int(1-X/x)E e5-d intLu+(d/2)int(1-X/x)Le',
        memory='Tin1=dMin5-(Deta Iin5+dJin5)/Hx;Tin2=(dU(X)-Deta)Min5-dSin5-Xd(Cpin5+axis5);no resets at U=0,modulation end or X=a',
        kernel='kappa=Xlog(X/a)-(X-a)<=ell^2/(2a);local quadratic swirl cost only,not zero total pressure or whole long interval smallness',
        inputs='A^-1<=a<=X<=A,values E,U,E^-1,Hx^-1<=A,C4 fields/moments/axis<=B=Ahi and lower pressure/shear C3<=B;FE,FU pure fifth/120;YE,YU y_eta4/24 remain actual inputs',
        remainder='truncated derivative moments C4<=4B,Pi<=2B,Pi_eta<=8B;Einv<=5A5B4,Hxinv<=5A6B4,Linv<=25/23<2;W21AB,Q156A7B5,N91AB2,R1<=2496A8B5,R2<=7280A7B6<=R0=10000A8B6',
        pressure='B1,B2 explicit positive equation4 linear in FE,FU and entering fifth bounds;|p1|4<=B+R0+10supB1,|p2|4<=B+R0+10A supB2;no actual input or gate evaluated',
        shear='highest angular fourth a,bs <=2A YE,YU +10A5B5;all lower denominator jets retained;not primitive a4,b4 or J5 inferred from J4',
        checks='twelve exact manufactured comparisons of both combined rows vs differentiated five genuine forward integrals,nonzero incoming/axis,pressure-drop controls,original-forward remainder and shear controls;not actual PDE samples',
        scope='same N widths roots support pressure gauge;native high inputs fourth scalar q5 moving-phase constants gate globalH23 fullpulse covariance means residual independentfoundations remainopen')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('volterra_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','volterra.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','fourth_pressure_input_boundary/input_bounds.py',
        'compact_pressure_bridge/bridge.py','higher_jet_budget/higher_bounds.py','annular_exterior_foundation/exterior_bounds.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Volterra record changed')
    if not all(expected['gates'].values()):raise ValueError('Volterra check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual native high inputs, compact gate and global H23 remain unverified.')
if __name__=='__main__':main()
