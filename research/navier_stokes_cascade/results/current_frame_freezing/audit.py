#!/usr/bin/env python3
"""Mixed current-frame freezing from value error; no retained-band or full pulse claim."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from mixed_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='mixed-current-frame-freezing-and-global-cone-cover-from-value-error'
DEPENDENCIES=('leading_frame_freezing','same_tilt_deficit','background_gap_certificate',
    'regional_weighted_gate','first_seminorm_bound','cutoff_cone_scale_gate','primary_covariance_foundation')
FLAGS=dict(current_endpoint_frame_comparison_uses_leading_jets_plus_two_value_errors=True,
    band_normalization_eight_and_current_orientation_omega_over_three_charged=True,
    mixed_mesh_and_background_freezing_terms_both_spend_cone_margin=True,
    global_compact_and_outer_sign_cover_retains_pointwise_full_pressure_sensitive_target=True,
    stress_magnitude_and_direction_derivatives_not_needed_for_this_value_comparison=True,
    full_corrected_width_pressure_viscosity_and_common_velocity_stress_Borel_inputs_retained=True,
    construction_choices_and_historical_evidence_unchanged=True,
    pure_mesh_uniform_current_frame_Lipschitz_bound_follows_from_C0_closeness=False,
    global_outer_relative_stress_smallness_or_magnitude_upper_bound_established=False,
    actual_cutoff_raw_deficit_endpoint_constants_and_retained_band_verified=False,
    arbitrary_generated_nonlinear_frames_satisfy_these_closeness_hypotheses=False,
    current_frame_derivative_and_all_actual_physical_pulse_gates_verified=False,
    complete_wave_covariance_axial_pressure_mean_and_residual_improvement_proved=False,
    independent_unconditional_core_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['supply_actual_unchanged_cutoff_norms_raw_deficit_and_endpoint_constants_and_verify_mixed_retained_scale_inequality',
    'verify_all_geometry_length_carrier_column_inverse_and_uniform_physical_pulse_value_derivative_support_gates',
    'establish_required_current_value_closeness_or_separate_bounds_for_each_generated_nonlinear_frame',
    'complete_full_wave_covariance_mean_corrections_including_axial_velocity_and_pressure_memory',
    'prove_complete_viscous_cutoff_interaction_residual_and_independently_check_retained_foundations']
def protocol():
    return dict(parent_commit='cd07a79f66c18f731a2a634eb35c7ee6ffffbd63',
        normalization='current band F,g divided by B=s^(-1-h);B^-1<8 on1/4<=s<=4;native endpoint error epsilon=8E;no cutoff or current perturbation differentiated',
        global_cover='compact epsilon<=rd gives orientation omega/3 and gammaB=omega^2*d0/(9kappa);outer current velocity equals leading and nonnegative angular stress increment improves orientation/cone without global relative-stress smallness',
        endpoint_floors='F>=f/2,-gtheta>=t/2,F<=Fbar+1,|g|<=Gbar+1;c2lower=d/[4(Fbar+1)(Gbar+1)],c2upper=2[(Gbar+1)^2+2(Fbar+1)(Gbar+1)]/(ft);rational cmin=min(1,c2lower),cmax=max(1,c2upper)',
        two_point='current differences bounded by leading JF*d+2epsilon,Jg*d+2epsilon;parent coordinate costs are linear separately,not their maximum',
        mixed_freezing='tau<=Cmesh ell^-6+8Cpert E;Cmesh=Dbox*max(LV(JF,Jg),LW(JF,Jg)),Cpert=max(LV(2,2),LW(2,2));use current endpoint floors and omega/3',
        precise_gate='8E<=rb;ell^6>=Dbox/dchart;Ract+max(16Lpsi E/omega^2,72kappa*(Cmesh ell^-6+8Cpert E)/omega^2)<=kappa;E=Cstar(M,N)*2^(-2h(ell-2))',
        later_witness='E<=min(rd/8,gammaB/(128Cpert));ell^6>=max(Dbox/dchart,16Cmesh/gammaB);original dyadic witness and exact sixth-root ceilings,not original first-band acceptance',
        pointwise_target='complete target stays at evaluation point in both coordinates;amplitude cancels;pressure-sensitive orientation and all compact/full-source premises retained;no representative stress frozen',
        obstruction='rational rotations+a,-a at distance a^2 have shear error<=8a and relative W error=4acw/[u(1-a^2)];value-only closeness does not bound slope uniformly;not actual PDE fields',
        diagnostics='192 exact current-two-endpoint controls with both target components,three tilts,eight signs,two distances/errors,two magnitudes;four oscillation controls,symbolic compact/outer identities,and auxiliary full-cost later witness',
        scope='conditional current-background VALUE freezing only;no actual scale,arbitrary generated frames,current derivatives,all physical pulses,complete covariance/residual,unconditional foundations or blowup certificate')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('current_freeze_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','mixed_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','leading_frame_freezing/freezing_bounds.py',
        'same_tilt_deficit/deficit_bounds.py','background_gap_certificate/gap_bounds.py',
        'cutoff_cone_scale_gate/gate_bounds.py','cutoff_shape_feasibility/shape_bounds.py',
        'actual_band_pulse_bounds/README.md','auxiliary_geometry_audit/README.md')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Current-frame freezing record changed')
    if not all(expected['gates'].values()):raise ValueError('Current-frame freezing check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual scale, arbitrary nonlinear frames and complete physical pulses remain unverified.')
if __name__=='__main__':main()
