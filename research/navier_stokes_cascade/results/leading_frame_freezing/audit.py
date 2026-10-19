#!/usr/bin/env python3
"""Explicit leading frame-only freezing; pointwise full target, no scale acceptance."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from freezing_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='explicit-leading-frame-freezing-from-global-jets-with-pointwise-full-target'
DEPENDENCIES=('same_tilt_deficit','global_swirl_envelope','global_axial_source_envelope',
    'leading_frame_foundation','source_envelope_edge_audit','primary_covariance_foundation')
FLAGS=dict(explicit_leading_reference_frame_freezing_constant_from_global_profile_jets=True,
    target_stress_is_kept_at_evaluation_point_not_frozen_at_representative=True,
    target_magnitude_cancels_without_direction_derivative_or_edge_amplitude_floor=True,
    pointwise_full_pressure_sensitive_orientation_and_cone_remain_required=True,
    enlarged_positive_chart_and_original_doubled_box_diameters_explicit=True,
    original_tilt_and_full_pressure_viscosity_deficit_gate_retained=True,
    construction_choices_and_all_historical_frozen_evidence_unchanged=True,
    actual_cutoff_raw_deficit_and_retained_first_band_verified=False,
    leading_freezing_constant_automatically_bounds_perturbed_or_nonlinear_frame=False,
    target_amplitude_can_be_frozen_with_uniform_relative_distance_error=False,
    full_pressure_source_vanishes_where_local_U_vanishes=False,
    unbounded_weighted_target_or_primary_root_quotients_made_bounded=False,
    column_length_inverse_carrier_and_all_actual_pulse_gates_verified=False,
    independent_unconditional_core_regional_or_global_cone_foundations_proved=False,
    complete_wave_covariance_mean_corrections_and_residual_improvement_proved=False,
    blowup_proved=False,formal_certificate=False)
REMAINING=['supply_actual_unchanged_cutoff_norms_raw_deficit_lower_bound_and_endpoint_profile_constants_and_check_retained_scale',
    'extend_frame_freezing_to_each_required_perturbed_or_nonlinear_frame_with_its_own_jets_orientation_and_value_floors',
    'verify_all_geometry_length_carrier_column_inverse_and_uniform_actual_pulse_value_derivative_support_gates',
    'complete_wave_covariance_mean_corrections_including_axial_velocity_and_full_pressure_memory',
    'complete_viscous_cutoff_and_all_interaction_residual_and_independently_verify_retained_foundations']
def protocol():
    return dict(parent_commit='b0f179b085aa0a64f281fd9ada47951c28e0e196',
        target='T(p)=m(p)n(p) stays at evaluation point in both V0,W0 and Vf,Wf; freeze only N,K,|c0|; both target components and complete pressure-sensitive orientation retained',
        frame_values='both endpoints f<=F<=Fbar,-gtheta>=t,|g|<=Gbar,cmin<=|c0|<=cmax,-N.n>=omega at evaluation point; native profiles after common positive band factor cancels',
        frame_jets='|DeltaF|<=JF*d,|Deltag|<=JG*d;LC=LD/(2ft)+Dbar*LP/(2f^2*t^2),LN=2JG/t',
        relative_coordinates='LV=cmax*LN/(cmin*omega)+LC/(2cmin^2);LW=cmax*(1+u)*LN/(u*omega); no stress or unit-direction derivative',
        positive_chart='T>=0;1/4<=s<=4,Xa<=X<=Xb at representative; enlarge to 1/8<=s<=8,Xa/2<=X<=2Xb;dchart=min(1/288,Xa/(4KX));KX=8b+288b^2,Keta=14',
        full_jets='global swirl native all-radius mixed envelope extends to 2Xb; axial derivatives extend by zero beyond actual Xv, not the pressure source; gz,X includes U_X/sqrt(2X)',
        boxes='retained halfwidth delta=ell^-6;original representative to k-fold box distance<2(k+1)delta;Dbox=4 original,6 doubled',
        freezing='Cfreeze=Dbox*max(LV,LW);ell^6>=max(Dbox/dchart,8Cfreeze/gamma0);frozen margin gamma0/2; other length,column,inverse and carrier gates separate',
        join='leading-only Cfreeze composes conditionally with same-tilt dreq and full Cstar; not a certificate for an arbitrary perturbed frame or original first band',
        countercontrol='m(y)=exp(-1/y^2),distance=y^3;amplitude ratio tends to exp(2),frame error zero; not actual PDE data',
        diagnostics='216 exact rational endpoint controls, three tilts, all sign choices, three distances, magnitudes1/10^40/10^200; symbolic identities and auxiliary complete global-jet/pressure/deficit join',
        scope='conditional leading-reference analytic estimate; no actual scale, all pulses, complete covariance/residual, unconditional foundations or blowup certification')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('leading_freeze_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','freezing_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','same_tilt_deficit/deficit_bounds.py',
        'global_swirl_envelope/swirl_bounds.py','global_axial_source_envelope/envelopes.py',
        'background_gap_certificate/gap_bounds.py','cutoff_cone_scale_gate/gate_bounds.py',
        'actual_band_pulse_bounds/README.md','auxiliary_geometry_audit/README.md')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Leading frame freezing record changed')
    if not all(expected['gates'].values()):raise ValueError('Leading frame freezing check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual scale and perturbed/nonlinear frame freezing remain unverified.')
if __name__=='__main__':main()
