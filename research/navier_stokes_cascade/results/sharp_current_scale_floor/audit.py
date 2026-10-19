#!/usr/bin/env python3
"""Sharp complete affine cutoff floor and mixed current-scale certificate; no actual acceptance."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from cost_floor import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='sharp-full-cutoff-cost-floor-and-coupled-current-scale-certificate'
DEPENDENCIES=('cutoff_shape_feasibility','cutoff_cone_scale_gate','current_frame_freezing',
    'same_tilt_deficit','first_seminorm_bound','leading_frame_freezing','regional_weighted_gate')
FLAGS=dict(sharp_infimum_of_all_four_fixed_affine_rows_derived=True,
    derivative_stress_crossing_and_stationary_branches_compared_exactly=True,
    equality_attainable_iff_constant_row_strictly_dominates_variable_floor=True,
    mixed_current_frame_budget_charges_both_mesh_and_background_freezing=True,
    shape_independent_required_deficit_floor_applies_to_this_sufficient_certificate=True,
    complete_pressure_shifted_viscosity_corrected_width_and_global_cone_cover_retained=True,
    historical_frozen_evidence_and_original_comparators_unchanged=True,
    full_infimum_is_attained_in_every_regime=False,
    actual_PDE_error_or_actual_raw_deficit_has_this_lower_bound=False,
    unchanged_cutoff_norms_and_historical_first_band_acceptance_verified=False,
    auxiliary_optimal_shape_replaces_old_cutoff_or_field=False,
    arbitrary_nonlinear_frames_or_all_actual_physical_pulse_gates_closed=False,
    complete_covariance_axial_pressure_means_and_residual_improvement_proved=False,
    independent_unconditional_core_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['bound_unchanged_cutoff_actual_norms_raw_deficit_and_endpoint_constants_then_check_full_mixed_retained_scale_gate',
    'verify_all_geometry_length_carrier_column_inverse_and_uniform_physical_scale_pulse_jets_supports',
    'establish_closeness_hypotheses_or_separate_bounds_for_every_required_nonlinear_frame',
    'complete_full_wave_covariance_and_axial_velocity_pressure_mean_corrections',
    'prove_complete_viscous_cutoff_interaction_residual_and_independently_check_retained_foundations']
def protocol():
    return dict(parent_commit='bbfe73f645de9184bc4cd4e52bd8cb8d64a099c2',
        fixed_program='Cstar=max(C0,D0+D1M,A0+A1M+A2N),C0=max(1,Rv),D1>0,A2>0,A1>=0;all other complete envelope inputs held fixed',
        shape_class='smooth K=1 before0,K=0 after1,0<=K<=1;capacity>1;existential restriction1<M<=2,N>M^2/(M-1),not reassigned norms of unchanged K',
        optimizer='P=A1+A2,Q=A2,S0=A0+A1+2A2,G0=D0+D1;ts=sqrt(Q/P);if g(ts)<=f(ts),V=S0+2sqrt(PQ);otherwise unique crossing t* in(0,ts),V=G0+D1t*',
        crossing='a=P-D1,b=G0-S0; t*=Q/b if a=0,else(b-sqrt(b^2-4aQ))/(2a)=2Q/(b+sqrt(b^2-4aQ))',
        attainment='I=max(C0,V);some smooth shape attains I iff C0>V;otherwise every smooth-shape cost strictly exceeds I',
        exact_decision='Cstar<=B possible iff B>I or(B=I and C0>V);exact rational quadratic-field sign arithmetic cross-checked with original capped quadratic',
        mixed_budget='H=gammaB/8-Cmesh/ell^6;Beff=S*min(rb/8,omega^2*d0/(16Lpsi),H/(8Cpert));require H>0,positive chart gate,and full shape decision at Beff',
        deficit_floor='dinf=max(16Lpsi*I/(omega^2*S),72kappa*Cmesh/(omega^2*ell^6)+576kappa*Cpert*I/(omega^2*S));strict for smooth shapes when C0<=V;certificate floor,not a physical error lower bound',
        decay_factor='exact mathematical S=2^(2h(ell-2));rational diagnostic S=2^floor(2h(ell-2)) is a conservative error-decay factor,not unrounded physical equivalence',
        pressure='full affine rows keep corrected tc,all future pressure moments,pointwise source,and shifted axial/swirl viscosities;complete pointwise target and compact/outer sign cover unchanged',
        diagnostics='10 exact optimizer regimes,414 independent exact budget comparisons,stationary/crossing/constant equality controls,mixed-bound obstruction and complete global-jet pressure composition;auxiliary fixtures only',
        scope='sharp smooth-shape result for fixed sufficient affine program;actual cutoff/deficit/first band,all pulses,nonlinear frames,full covariance/residual and unconditional foundations remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('sharp_current_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','cost_floor.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','cutoff_shape_feasibility/shape_bounds.py',
        'cutoff_cone_scale_gate/gate_bounds.py','current_frame_freezing/mixed_bounds.py','same_tilt_deficit/deficit_bounds.py',
        'leading_frame_freezing/freezing_bounds.py','first_seminorm_bound/seminorms.py','actual_band_pulse_bounds/README.md')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Sharp complete cutoff/current-scale record changed')
    if not all(expected['gates'].values()):raise ValueError('Sharp complete cutoff/current-scale check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual cutoff and retained-scale acceptance remain unverified.')
if __name__=='__main__':main()
