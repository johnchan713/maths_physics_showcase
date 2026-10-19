#!/usr/bin/env python3
"""Same heat-root third jets and full leading swirl input; no scale acceptance."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from swirl_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='same-heat-root-third-jets-and-global-swirl-input-envelope'
DEPENDENCIES=('global_axial_source_envelope','global_axial_angular_bound','root_foundation',
    'annular_exterior_foundation','higher_jet_budget','cutoff_repair_derivative_budget','first_seminorm_bound')
FLAGS=dict(heat_target_C3_derived_from_positive_Gamma_integral_not_relabelled_from_C1=True,
    full_infinite_heat_angular_moment_h_times_one_over_h_retained=True,
    same_heat_root_identified_in_common_pointwise_uniqueness_ball=True,
    full_reference_swirl_J3_includes_heat_compensation_angular_edits_and_terminal_wait=True,
    exact_original_normalized_two_unit_terminal_taper_retained=True,
    Gamma_jets_use_only_nonnegative_argument_and_one_sided_pole_limits=True,
    global_F_radial_and_shifted_swirl_inputs_join_axis_compact_and_later_reference=True,
    all_primitive_leading_first_stress_inputs_now_explicit_relative_to_retained_lemmas=True,
    original_roots_profiles_supports_cutoffs_parameters_gauges_and_both_mean_memories_unchanged=True,
    historical_frozen_evidence_preserved=True,
    new_larger_C3_root_ball_is_inside_old_C1_ball=False,
    actual_cutoff_shape_L1_L2_and_full_Cstar_tilt_scale_acceptance_verified=False,
    actual_huge_profile_or_physical_pulses_numerically_resolved=False,
    retained_core_and_regional_premises_independently_proved_unconditional=False,
    unconditional_pulse_complete_residual_improvement_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['identify_or_bound_actual_unchanged_coefficient_cutoff_shape_L1_L2_then_verify_full_same_field_Cstar_tilt_scale_and_other_acceptance_gates',
    'prove_uniform_actual_physical_pulse_value_derivative_inverse_gap_and_support_estimates_on_every_scale',
    'retain_complete_wave_covariance_axial_velocity_pressure_and_mean_terms_in_full_viscous_interaction_residual',
    'independently_verify_retained_exact_core_and_regional_foundations_before_unconditional_claims']
def protocol():
    return dict(parent_commit='ca0e6a3aec5b491bfd0069418261a092d8fdd03e',
        norms='factorial C_eta^3 and mixed total J3, including one-sided eta endpoint jets',
        heat_target='|Hcal-1|3<=512h/X; complete raw heat rows <=(3072h eK^2/XK,3072h eK^2,1024sqrt(2)eK sqrt(XK)); normalized target <=8192h^4/Xstar',
        same_heat_root='B inverse<15,Q<29; root C3<=245760h^4/Xstar<2^18h^4/Xstar; same old root in common pointwise ball, not C3 subset of old C1 ball',
        original_terminal='fo(y)/fo(0),fo=1-(h/1000)*sigma((3-y)/2),fo(0)=1-h/1000; same width two and orientation',
        reference='scalar J3<8192Ahi, angular shape J3<2^34; four actual factors J3<2; later native E J3<2^60Ahi<H^2',
        global_F='native E J3<=H^70 after axis; MF=2(1+Xa^-1)H^70; F3=max(Kaxis,MF),FX2=max(Kaxis,MF/Xa),FXX1=max(Kaxis,3MF/Xa^2)',
        shifted_swirl='CF=64Xc^2 FXX1+704Xc FX2+1664F3, all coordinate derivatives retained; compose same first source and compact Cstar interface',
        diagnostics='ten positive Gamma integrals at 45 digits, displayed point jets rounded to 20 significant digits; not actual-profile samples',
        scope='conditional leading-field input completion; actual coefficient cutoff shape and full Cstar/tilt/scale gates, physical pulse/covariance/residual and blowup remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('global_swirl_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','swirl_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','global_axial_source_envelope/envelopes.py',
        'root_foundation/root_bounds.py','outer_pressure_pilot/schedule.py','heat_exterior_audit/construction.py','first_seminorm_bound/seminorms.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Global swirl/heat record changed')
    if not all(expected['gates'].values()):raise ValueError('Global swirl/heat check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual coefficient cutoff shape and full same-field stress acceptance remain unverified.')
if __name__=='__main__':main()
