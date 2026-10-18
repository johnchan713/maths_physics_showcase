#!/usr/bin/env python3
"""Explicit background value gap, not independent verification of every leading input."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from gap_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='explicit-background-value-gap-certificate-relative-to-leading-and-weighted-jets'
DEPENDENCIES=('leading_frame_foundation','weighted_background_edges','weighted_stress_foundation',
 'physical_pulse_foundation','canonical_primary_foundation','sourced_inverse_foundation','current_residual_ledger_audit')
FLAGS=dict(denominator_free_growth_discriminant_and_exact_c_squared_retained=True,
 explicit_rational_radius_replaces_unspecified_background_continuity_constants=True,
 same_existing_positive_tilt_preserved_with_explicit_weaker_cone_margin=True,
 actual_weighted_relative_stress_input_required_at_zero_edges=True,
 value_domain_selected_once_not_by_high_derivative_order=True,
 full_four_Q_chart_cost_and_every_previous_tail_gate_retained=True,
 canonical_seed_route_already_completed_not_counted_as_new_result=True,
 no_frozen_phase_frame_carrier_root_width_or_Borel_schedule_redefined=True,
 full_Pi1_Pi2_pressure_A0_A1_and_exact_residual_ledger_preserved=True,
 exact_finite_controls_not_actual_giant_profile_values_or_continuum_proof=True,
 recursive_byte_verification_separate_from_fresh_dependency_replays=True,
 actual_Cval_and_huge_leading_constants_numerically_evaluated=False,
 all_primitive_core_regional_leading_weighted_inputs_independently_proved_here=False,
 historical_first_band_or_existing_qstar_passes_new_sufficient_gate=False,
 actual_background_shear_frame_substituted_for_frozen_pulse_frame=False,
 generic_absolute_stress_error_suffices_at_flat_edges=False,
 failure_of_sufficient_test_proves_actual_cone_failure=False,
 all_fixed_jets_small_on_one_domain_or_uniform_in_order=False,
 unconditional_pulse_residual_single_field_force_or_blowup_proved=False,formal_certificate=False)
REMAINING=['independently_discharge_actual_regional_small_root_value_and_weighted_error_inputs_in_the_explicit_Cval_q_power_gate',
 'retain_historical_first_band_and_old_force_limits_and_join_any_later_distinct_realization_with_its_own_controlled_schedules',
 'do_not_repeat_completed_canonical_seed_route_or_promote_a_background_value_certificate_to_an_unconditional_complete_residual_or_blowup_theorem']

def protocol():
    return dict(parent_commit='165b75afee625f8ec597e117e4986452907e7a3c',
      inputs='f<=F0<=Fbar, minus_gtheta>=t, norm_g0<=Gbar,D0>=d; H0>=omega; same u cone ratio<=1-gamma; normalized error<=Cval*q^(2h)',
      discriminant='D=norm_g^2+2F*gtheta; lambda2=-2F*gtheta*D/norm_g^2; c2=D/(-2F*gtheta)',
      radius='min(1/2,f/2,t/2,d/(2LD),omega/(2LH),mphi/(2Lphi)) with explicit LD,LP,Lc,LH,Lphi',
      cone='Phi=u^2 H^2-(1+u^2)c^2 J^2; actual margin omega^2*gamma*(2-gamma)/16, existing u unchanged',
      floors='F>=f/2,minus_gtheta>=t/2,D>=d/2,lambda2>=f*t*d/[4(Gbar+1)^2], c2>=d/[4(Fbar+1)(Gbar+1)]',
      domain='Cval*q^(2h)<=radius; whole dyadic chart needs ell>=2+log2(Cval/radius)/(2h), combined with every previous gate',
      quantifiers='one value domain, every fixed smooth jet with order dependent constants; real smoothness not complex heat/modulation analyticity',
      memory='full Pi1/Pi2 pressure, axial velocities, A0/A1, normal wave pressure, temporal Gaussian vs regular spatial cutoff distinction all unchanged',
      scope='conditional explicit background value-gap interface only; original first band, actual constants, primitive leading/core premises and unconditional theorem not certified')

@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('background_gap_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','gap_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py']
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
    e=generate()
    if json.dumps(e,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Background-gap record changed')
    if not all(e['gates'].values()):raise ValueError('Background-gap gate failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates; explicit real-domain spectral and same-tilt cone bounds")
    print('Leading/weighted inputs, actual first-band gate and unconditional theorem remain unpromoted.')
if __name__=='__main__':main()
