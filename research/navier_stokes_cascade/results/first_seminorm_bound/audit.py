#!/usr/bin/env python3
"""Quantitative first-coefficient input bound, not actual primitive envelope evaluation."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from seminorms import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='collar-linear-repair-full-pressure-first-seminorm-bound'
DEPENDENCIES=('first_correction_construction','annular_first_correction','coefficient_support_foundation',
              'weighted_stress_foundation','regional_weighted_gate','background_gap_certificate')
FLAGS=dict(collar_Picard_and_angular_Cauchy_bounds_propagated_through_actual_linear_repair=True,
 radial_first_two_jets_use_regular_system_not_radial_analyticity=True,
 full_global_Omega0_over_X_moment_memory_and_nonconstant_amplitude_normalizers_retained=True,
 cutoff_and_logarithmic_bump_jets_and_E_to_F_inverse_radius_costs_charged=True,
 complete_first_pressure_value_and_axial_derivative_and_shifted_leading_axial_viscosity_kept=True,
 exact_inner_zero_stress_avoids_dividing_at_flat_Xa=True,
 explicit_finite_Cstar_composes_with_previous_regional_gate=True,
 previous_outer_heat_sign_and_seed_results_not_recounted=True,
 actual_Ccol_tau_primitive_leading_bounds_and_combined_gate_evaluated=False,
 core_only_C16_operator_bound_valid_on_actual_activated_collar=False,
 Picard_majorant_bounded_by_a_fixed_polynomial_in_C=False,
 first_stress_bound_is_complete_wave_covariance_or_complete_residual_bound=False,
 all_other_root_value_stress_magnitude_and_residual_gates_discharged=False,
 historical_first_band_old_field_or_old_force_recertified=False,
 unconditional_pulse_complete_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['justify_quantitative_actual_collar_jet_global_pressure_source_and_bump_cutoff_envelopes_and_compact_leading_stress_gap_inputs',
 'verify_2h_times_ell_minus_two_at_least_log2_Cstar_over_regional_radius_together_with_all_other_accepted_tail_gates',
 'retain_full_pressure_and_all_scale_pulse_covariance_mean_and_residual_obligations_and_do_not_recertify_historical_first_band']

def protocol():
    return dict(parent_commit='b29ca142b37c4200d5a90b824cb2824f70cc0688',
       collar='M=(z+beta)exp(beta); z=2*Ccol*a_col; beta=e*z^2/Delta; Delta+d<sigma; actual collar Ccol NOT core C^16',
       jets='J=M(1+1/d); Jeta=M(1/d+2/d^2); H=Ccol(1+1/d)(J+Jeta+1); Jx=J/(2a); Jxx=H/(4a^2)+J/a^3; Xa=a^2',
       repair='five C_eta^1 ordinary moments include int_0^Xv |Omega0/X|_1 and int_0^Xv X|Omega0/X|_1; inverse 800; reciprocal jet cost 4/e_star; lambda^-1 kept',
       radii='Xa=a^2; Xc=c^2; sqrt(2X) bounded below by a and above by 2c; rho supplied as rho_root^2; p lower patch radius and rp^2<=2p',
       pressure='Pi1(0)=0; Pi1_X=2F0F1-Omega0/(2X); complete pressure Z_-1 retained; no direct pressure cutoff or zero radial-memory assumption',
       stress='T1=0 through Xinner>Xa; Ttheta=-int x*a_theta dx/X; Tz=-int a_z dx/sqrt(2X); tau=min_[Xinner,Xc] |T0|>0',
       gate='Cstar from derived first seminorm bounds and previous shared Borel budget; sufficient whole-chart inequality 2h(ell-2)>=log2(Cstar/r_Kc)',
       scope='analytical finite input bound conditional on primitive leading/collar lemmas; diagnostic constants are not selected-profile values or historical first-band acceptance')

@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('first_seminorm_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','seminorms.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py',HERE.parent/'first_correction_construction/system.py',
        HERE.parent/'annular_first_correction/extension.py',HERE.parent/'coefficient_support_foundation/support_bounds.py',
        HERE.parent/'regional_weighted_gate/weighted_gate.py',HERE.parent/'physical_residual_budget/residual.py',
        HERE.parent/'physical_residual_budget/review.py']
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('First seminorm record changed')
    if not all(expected['gates'].values()):raise ValueError('First seminorm check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual primitive envelopes and combined gate remain unverified.')
if __name__=='__main__':main()
