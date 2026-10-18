#!/usr/bin/env python3
"""Reproduce signed full-pressure source moments; scope remains conditional."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from moments import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='full-support-signed-pressure-source-moments-with-axis-viscosity'
DEPENDENCIES=('collar_operator_envelope','first_seminorm_bound','coefficient_support_foundation',
    'global_pressure_moment_foundation','annular_first_correction')
FLAGS=dict(full_nonconservative_radial_source_retained=True,
    exact_conservative_radial_flux_and_two_signed_moment_identities_derived=True,
    unweighted_axis_viscosity_term_retained_and_weighted_total_cancellation_distinguished=True,
    full_future_axial_support_and_functional_zero_M_premise_retained=True,
    finite_endpoint_corrections_prevent_truncated_future_memory=True,
    signed_total_norms_replace_stronger_absolute_integrals_only_for_repair=True,
    exact_auxiliary_control_rejects_promoting_signed_bound_to_absolute_integral=True,
    repair_source_bound_uses_global_zero_radial_third_angular_axial_norm=True,
    same_linear_inverse_amplitude_angular_normalizers_and_pressure_gauge_preserved=True,
    full_I4_Pi1_Pi2_pressure_memories_not_erased=True,
    historical_frozen_evidence_preserved=True,
    actual_global_axial_angular_norm_quantified_from_all_original_stages=False,
    signed_moment_bound_is_absolute_integral_or_pointwise_source_bound=False,
    full_global_first_stress_coefficient_inputs_discharged=False,
    actual_unchanged_tilt_and_combined_scale_gate_verified=False,
    historical_first_band_or_original_full_field_force_recertified=False,
    unconditional_actual_pulses_full_residual_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['quantify_actual_global_KU_through_all_later_leading_stages_not_just_activation_collar',
    'quantify_separate_pointwise_w0_global_shifted_axial_sources_and_cutoff_repair_bump_envelopes_for_Cstar',
    'verify_same_unchanged_tilt_and_all_scale_freezing_geometry_inverse_acceptance_gates',
    'retain_complete_covariance_axial_velocity_and_both_mean_pressure_memories_in_full_residual_closure']

def protocol():
    return dict(parent_commit='a70349e4b48b41c08033921824951ce286766195',
        source='same retained regular six-variable first-correction source, not final manuscript theorem',
        domain='full leading axial support [0,Xv] and eta in [-1,1]; regular v=V/X and exact functional M_infinity=0',
        flux='J=(1-2eta U)V/L+X*v^2-2nu V_X; w0=J_X+(Deta v_eta+d(Uv)_eta-2Aeta Uv)/L+v^2/2',
        I0='2nu v(0)+(Deta a0_eta+d b0_eta-2Aeta b0)/L+c0/2',
        I1='(-a1+Deta a1_eta+d b1_eta+2Deta b1)/L-c1/2',
        integrals='a_j=int X^j v,b_j=int X^j Uv,c_j=int X^j v^2, through full Xv',
        norm='C_eta^k=sum sup|partial_eta^j|/j!; Q0,Q1 given by README (8),(9), not int absolute w0',
        axial_input='KU=sup_[0,Xv] |U|_3; u2<=KU,v2<=36KU; actual full-stage KU not evaluated',
        support='arbitrary finite endpoint requires both explicit boundary corrections; axis term still retained',
        diagnostics='80 exact manufactured polynomial cases, 320 scalar moment/first-angular-jet identities; 18 floating retained-source comparisons at tolerance 2e-10',
        scope='conditional source-moment simplification and repair input only; no complete first stress, actual pulse, full residual or blowup certificate')

@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('signed_pressure_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','moments.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py',HERE.parent/'first_correction_construction/system.py',
        HERE.parent/'first_seminorm_bound/seminorms.py',HERE.parent/'physical_residual_budget/residual.py',
        HERE.parent/'coefficient_support_foundation/conservative_checks.py']
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Signed pressure source record changed')
    if not all(expected['gates'].values()):raise ValueError('Signed pressure source check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual global angular envelope and complete residual remain open.')
if __name__=='__main__':main()
