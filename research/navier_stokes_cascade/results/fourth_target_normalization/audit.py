#!/usr/bin/env python3
"""Reproduce tame fourth-angular normalized five-row repair target bounds."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from target_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='tame-fourth-angular-density-moments-and-normalized-target-primitive-gate'
DEPENDENCIES=('fourth_angular_axial_interface','higher_jet_budget','modulation_repair_all_order',
    'stress_realization_audit','coefficient_support_audit')
FLAGS=dict(tame_fourth_exponential_and_all_five_complete_density_differences_bounded=True,
    actual_normalized_target_Z4_bound_affine_in_two_highest_primitive_jets_derived=True,
    distinct_five_row_normalizers_fourth_derivatives_inverse_lambda_and_geometry_retained=True,
    same_original_preconditioner_checks_manufactured_fourth_target_jets=True,
    actual_primitive_B4_and_target_Z4_eliminated_into_one_sufficient_majorant_gate=True,
    moving_inverse_phase_fourth_derivatives_not_frozen_or_inferred_from_C3=True,
    full_incoming_moments_and_pressure_axis_datum_cancel_only_in_same_profile_differences=True,
    roots_width_frequency_supports_historical_evidence_and_future_mean_pressure_memory_preserved=True,
    actual_moving_phase_fourth_primitive_jets_a4_b4_and_substituted_gate_evaluated=False,
    actual_global_leading_U_C4_first_U_C3_or_second_H23_closed=False,
    higher_pressure_w1_viscosities_or_P2_lower_F1_square_replaced_by_leading_moment_densities=False,
    retained_scale_acceptance_and_complete_physical_pulse_derivatives_supports_verified=False,
    full_covariance_axial_pressure_means_and_viscous_cutoff_interaction_residual_proved=False,
    independent_unconditional_core_Picard_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['quantify_actual_two_fourth_slow_primitive_jets_a4_b4_including_variance_root_moving_inverse_phase_and_periodic_integration',
    'verify_substituted_compact_C4_majorant_gate_at_same_retained_parameters_and_derive_actual_global_leading_U_C4',
    'supply_finalized_first_U_C3_actual_coefficient_cutoff_jets_and_complete_second_pressure_repair_for_H23',
    'prove_retained_scale_acceptance_all_physical_pulse_jets_supports_complete_wave_covariance_axial_pressure_means_and_full_viscous_cutoff_interaction_residual',
    'independently_check_retained_core_sparse_Picard_and_regional_foundations']
def protocol():
    return dict(parent_commit='4f093b36c10abea2fb4c8e6b31bb522418df7963',
        inputs='same nativeE,U Ceta4<=Ahi,compact A^-1<X<A,length<A;actual primitivejointJ3<=H,N>=H,selectedN>H^32 unchanged;actuala4,b4 are fourth normalized eta jets at fixed external phase with full internal moving inverse phase',
        tame='E_N=E exp(a),U_N=U+b;a=Acal/N,b=Bcal/N;lowerC3<=t=H/N<=1;exp(a)-1 C4<=(12H+3a4)/N,exp(2a)-1 C4<=(114H+18a4)/N;b^2 C4<=(4H+2b4)/N;b exp(a) C4<=(26H+3a4+3b4)/N',
        densities='N times full C4 differences:DU=H+b4,DE=Ahi(12H+3a4),DU2=Ahi(6H+4b4),DE2=Ahi^2(114H+18a4),DUE=Ahi^2(38H+6a4+3b4);both squares and complete product retained;linear in highest jets',
        raw='mM=A DU,mJ=2A^2 DUE,mI=2A^2 DE,mS=A(DU2+DE2/2),mCp=A^2 DE2/2;complete forward incoming differences;Pi_axis unchanged,deltaPi=deltaCp notPi_axis=0',
        normalizers='samee=e_star/(1+eta^2),rhop>0,srho>=1/(sqrt2 rhop^1.5),0<lambda<=.01;wM=5/(rhop e_star),wJ=25srho/e_star^2,wI=5srho/e_star,wS=25/(rhop e_star^2),wCp=25/e_star^2;all C4 derivatives charged',
        transform='normalize FIRST then(m,(m-j)/lambda,i,-s,c);row C4 boundmax[wM mM,(wM mM+wJ mJ)/lambda,wI mI,wS mS,wCp mCp];componentwise coefficientmax givesZH,Za,Zb withZ4<=ZH H+Za a4+Zb b4',
        gate='parentB4<=H+b4,r=2000H^16/N,Kknown=Ahi+6400Ar+960000000000Ar^2;(3200000A Za)a4+(1+3200000A Zb)b4<=N(K-Kknown)-(1+3200000A ZH)H;equivalent only to substituted sufficient majorant,not necessity for actual field',
        diagnostics='fifteen manufactured density/target fourth jets ateta endpoints/interior and3phases;direct differentiation vsLeibniz andoriginal precondition_target;moving lift fourth9/512;symbolicexponentialpartitions;not actual loop/PDE or interval certification',
        pressure='leading modulation five-density target ONLY;not replacement forfullhigherw1 shifted viscosities,P2 F1^2 orPi1/Pi2 I4 memory;all five exact repair rows restore eta-functions andsame future pressure',
        scope='actuala4,b4,substitutedgate,globalleadingC4,firstC3,H23,cutoff/scale,pulses,covariance means,residual andindependentfoundations remainopen')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('fourth_target_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','target_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','fourth_angular_axial_interface/angular_bounds.py',
        'stress_realization_audit/repair.py','modulation_repair_all_order/all_order_maps.py','higher_jet_budget/higher_bounds.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Fourth target record changed')
    if not all(expected['gates'].values()):raise ValueError('Fourth target check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual moving-phase a4/b4, compact gate and global H23 remain unverified.')
if __name__=='__main__':main()
