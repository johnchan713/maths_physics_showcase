#!/usr/bin/env python3
"""Reproduce the fourth-angular reference bound and explicit compact gate."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from angular_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='same-reference-root-fourth-angular-axial-bound-and-explicit-compact-input-gate'
DEPENDENCIES=('axis_second_coefficient','global_axial_angular_bound','axial_root_foundation',
    'modulation_repair_all_order','higher_pressure_source_moments','higher_jet_budget',
    'cutoff_higher_jet_interface','collar_operator_envelope')
FLAGS=dict(same_reference_angular_root_C4_identified_with_unchanged_accepted_root=True,
    complete_power_reference_energy_C4_includes_angular_edits_wait_infinite_tail=True,
    same_strength_fourth_normalized_jets_retain_affine_end_corrections_and_odd_cross_terms=True,
    later_reference_axial_pulse_C4_bound110_derived_relative_to_retained_roots=True,
    compact_fourth_root_jet_uses_actual_transformed_target_and_pointwise_inverse=True,
    explicit_compact_B4_Z4_sufficient_majorant_gate_derived_without_frequency_reselection=True,
    conditional_global_fourth_axial_cover_and_full_signed_pressure_interface_supplied=True,
    historical_roots_width_frequency_supports_full_pressure_and_mean_memories_preserved=True,
    actual_compact_loop_B4_normalized_target_Z4_and_gate_evaluated=False,
    finalized_actual_global_leading_U_C4_or_first_U_C3_quantified=False,
    finalized_second_H23_and_unchanged_coefficient_cutoff_scale_acceptance_closed=False,
    heat_modified_R_alone_identified_with_power_reference_R=False,
    full_physical_perturbation_pulse_derivatives_supports_covariance_means_and_residual_closed=False,
    independent_unconditional_core_Picard_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['quantify_actual_compact_fourth_angular_primitive_B4_and_fully_normalized_five_row_target_Z4_including_moving_phase_and_inverse_lambda_costs',
    'verify_compact_majorant_gate_and_global_leading_U_C4_without_reselecting_frequency_or_root',
    'supply_finalized_global_first_U_C3_actual_coefficient_cutoff_jets_and_complete_second_repair_for_H23',
    'prove_retained_scale_acceptance_all_physical_pulse_derivatives_supports_full_wave_covariance_axial_pressure_means_and_complete_viscous_cutoff_interaction_residual',
    'independently_check_retained_core_sparse_Picard_and_regional_foundations']
def protocol():
    return dict(parent_commit='5ad730f469503d4d817e3abcd6d58390d8281152',
        norms='factorial angular Ck onclosed[-1,1];f=1/(1+eta^2),a=log(2/(1+eta^2));pole distances>=1 givef_C4<=5,a_C4<31/6<8',
        angular_root='same exact discrepancy formula impliesdeltaI_C4<2^21,z_C4<=2^29lambda^30,cang_C4<=2^32lambda^30;same pointwise uniqueness ball,no reselected root',
        energy='g_theta fourth derivative has all three terms,raw caps(1,4,12,56,408),C4<64;unedited1024T+3200T=4224T;p_squared_C4=25;signed edits30c+1000c^2<1;complete power R_C4<=8192T',
        heat='R is power reference;heat compensation and later heat edit have zero netS change as eta/trial-strength functions;do not identifyheat-modifiedR alone withR',
        strength='same k Amp^2+B Amp+D with both affine end bumps;|B|4<=160lambda^41,|D+C|4<=epsilon=8192Tlambda<2^-20;inverseJac<3;normalizedu1..u4<=12epsilon;Amp_C4<2,ci_C4<=4lambda^20;no parity assumed',
        later_pulse='f_C4<=5,R0<=11,main andend supportsdisjoint,beta<=64,eb exponential<=1;reference axial pulseC4<=110,taper<=8,laterU0=0',
        compact='same originalU_C4<=Ahi,U_N=U+Bcal(y,eta,Ny)/N;actual B4 full angular norm;actual Z4>=N sup||z_N^(4)||inf/4! includes all incoming moments,e-normalizers,inverse-lambda row andtheir derivatives',
        same_repair='r=2000H^16/N boundsd0..d3;constantB,Q,JacInverse<=2000;d4<=2000Z4/N+6e8r^2;sumc_C4<=4r+d4;disjointunitmasswidth.1bumps<=320,e_C4<=5A',
        gate='Kcompact4<=Ahi+B4/N+1600A(4r+2000Z4/N+6e8r^2);majorant<=K iff B4+3200000A Z4<=N(K-Ahi-6400Ar-960000000000Ar^2);not necessary for truefield or new schedule',
        global_pressure='Kaxis4=K(ell_axis)sumdelta^-j j0..4;Kleading4=max(Kaxis4,Kcompact4,110);supply actual globalKfirst3 then regular_v_bound givesv0_C2,v1_C2,v0_C3;full k1 signedI0,I1 through originalXv retainboth shifted viscosities andpressurememory',
        diagnostics='five manufactured positive quadratic roots all four exact normalized recurrences;abstractn^-3cos(neta)sin(phi) C3<=8/3 butnormalizedfourth n/24 isnot actual loop/PDE;tests not continuum proof',
        scope='actualB4,Z4,globalleadingC4,globalfirstC3,cutoff/finalH23,scale,physicalpulses,covariance means,residual andindependentfoundations remainopen')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('fourth_angular_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','angular_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','global_axial_angular_bound/angular_bounds.py',
        'collar_operator_envelope/envelope.py','higher_pressure_source_moments/moments.py',
        'modulation_repair_all_order/all_order_maps.py','cutoff_higher_jet_interface/jet_bounds.py','axial_root_foundation/axial_bounds.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Fourth-angular record changed')
    if not all(expected['gates'].values()):raise ValueError('Fourth-angular check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual B4/Z4, global U jets and H23 remain unverified.')
if __name__=='__main__':main()
