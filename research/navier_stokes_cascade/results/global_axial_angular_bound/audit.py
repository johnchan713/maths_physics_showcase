#!/usr/bin/env python3
"""Global axial third-angular input relative to retained exact regional lemmas."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from angular_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='global-third-angular-axial-envelope-and-signed-pressure-source-input'
DEPENDENCIES=('signed_pressure_source_moments','collar_operator_envelope','higher_jet_budget',
    'axial_root_foundation','root_foundation','matching_root_foundation')
FLAGS=dict(angular_target_C3_derived_not_relabelled_from_C1=True,
    same_angular_root_identified_in_common_pointwise_uniqueness_ball=True,
    complete_post_pulse_energy_C3_includes_all_edits_wait_and_infinite_tail=True,
    full_affine_end_energy_and_non_even_strength_terms_retained=True,
    normalized_third_strength_jets_and_actual_composed_end_coefficients_bounded=True,
    axis_compact_modulation_repair_and_reference_pulse_domains_joined=True,
    explicit_global_KU_and_two_signed_pressure_source_bounds_supplied_conditionally=True,
    original_width_frequency_roots_supports_pressure_and_both_mean_memories_unchanged=True,
    historical_frozen_evidence_preserved=True,
    actual_giant_parameter_values_numerically_materialized=False,
    retained_core_and_compact_estimates_independently_proved_unconditional=False,
    signed_source_bound_is_absolute_integral_or_pointwise_residual_bound=False,
    complete_first_stress_cutoff_bump_and_shifted_axial_inputs_closed=False,
    unchanged_tilt_combined_scale_or_historical_first_band_verified=False,
    unconditional_actual_pulses_full_residual_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['quantify_separate_pointwise_w0_shifted_axial_sources_leading_derivatives_and_actual_cutoff_repair_bump_envelopes_in_Cstar',
    'verify_same_unchanged_tilt_and_all_scale_freezing_geometry_inverse_acceptance_gates',
    'retain_complete_wave_covariance_axial_velocity_and_both_mean_pressure_memories_in_full_residual_closure',
    'independently_discharge_remaining_retained_core_regional_and_actual_pulse_foundation_premises']

def protocol():
    return dict(parent_commit='0ba0d99a92c7e6de0158f91a11a49709f805ec64',
        norm='factorial C_eta^3 on closed [-1,1], one-sided endpoint jets; no derivative-order-uniform claim',
        parameters='same T>=64,lambda=exp(-4T),C>1e9,t=C^-20,Ahi=(C^4096)^256 and unchanged N',
        angular_root='|deltaI|3<2^21, |z|3<=2^29 lambda^30, |c_ang|3<=2^32 lambda^30; same pointwise root',
        complete_energy='|R|3<=8192T, including angular energy increments, infinite reference tail and unchanged waiting; later heat compensation has exact net S zero',
        strength='epsilon=8192Tlambda<2^-20; normalized jets <=(9,279,15939)epsilon, |Amp|3<2, |c_end|3<=4lambda^20',
        pulse='|f|3<16,R0<=11,disjoint supports imply |U|3<=352; no actual giant profile sampled',
        axis='collar envelope proof at ell=tc/2<t/20; Kaxis=K(ell)*(1+delta^-1+delta^-2+delta^-3)',
        global_bound='KU<=max(Kaxis,2Ahi,352); use same Xv and u2=KU,v2=36KU in signed Q0,Q1',
        scope='conditional explicit global axial/source-moment input; full first stress, same tilt, scale gates, actual pulses and complete residual remain open')

@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('global_axial_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','angular_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py',HERE.parent/'outer_pressure_pilot/schedule.py',
        HERE.parent/'pulse_moment_audit/pulse.py',HERE.parent/'signed_pressure_source_moments/moments.py',
        HERE.parent/'collar_operator_envelope/envelope.py']
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Global axial envelope record changed')
    if not all(expected['gates'].values()):raise ValueError('Global axial envelope check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; full first-stress, tilt, scale and residual acceptance remain open.')
if __name__=='__main__':main()
