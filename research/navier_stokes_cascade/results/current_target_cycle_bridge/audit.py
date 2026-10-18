#!/usr/bin/env python3
"""Conditional exact current-target recurrence, not an unconditional residual proof."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from target_cycle import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='conditional-current-target-cancellation-cycle-bridge'
DEPENDENCIES=('current_target_moment_audit','source_envelope_edge_audit','canonical_finite_cycle_join',
 'canonical_initial_mean_join','compact_mean_foundation','complex_mean_relaxation','complex_signed_covariance_update')
FLAGS=dict(current_target_estimate_derived_by_exact_state_difference_not_raw_field_bound=True,
    generic_signed_target_cancellation_keeps_all_six_covariances_and_changed_pressure=True,
    moving_compensator_derivative_commutator_retained=True,
    actual_zero_Haar_fast_curl_has_exactly_zero_linear_Haar_response=True,
    fast_old_new_correlations_and_recomputed_postfast_defects_retained=True,
    slow_current_pressure_time_viscosity_and_old_Haar_means_retained=True,
    all_half_line_target_recurrence_relative_to_reviewed_current_map_inputs=True,
    accepted_canonical_initial_target_routes_to_fresh_target_seed=True,
    exact_linear_nonzero_angular_tail_routing_keeps_full_cut_velocity_covariance=True,
    same_complete_base_Pi1_Pi2_actual_A0_A1_fixed_width_and_corrected_domain_retained=True,
    historical_sources_evidence_and_parent_smooth_countercontrol_preserved=True,
    manufactured_tensor_updates_not_actual_positive_wave_realizations=True,
    initial_canonical_wave_mean_and_inverse_estimates_independently_reproved_here=False,
    unconditional_full_residual_and_all_analytical_premises_closed=False,
    uniform_in_stage_jet_width_or_lambda_constants_proved=False,
    arbitrary_flat_tails_forward_inverted=False,
    actual_giant_fields_first_band_or_old_force_acceptance_evaluated=False,
    unconditional_pulse_or_blowup_proved=False,formal_certificate=False)
REMAINING=['audit_complete_nonzero_angular_residual_increment_and_exact_pressure_cutoff_tail_reclassification_with_every_current_interaction_and_physical_jet_support_cost',
 'independently_discharge_retained_actual_canonical_initial_wave_mean_background_core_and_inverse_map_estimates_before_unconditional_claims',
 'do_not_identify_the_accepted_distinct_tail_witness_with_historical_first_band_fields_or_forces_or_use_an_uncontrolled_stage_diagonal']

def protocol():
    return dict(parent_commit='f370f9155a7286d613ad4e00ba5e22a3ddf43b0d',
        target='Xe=-Ie Qe bar(Ee); dXe=-Ie Qe bar(dEe), all CURRENT conservative equations and both velocity moments',
        wave='dXtheta=-Drtheta-I2 Q2 Dz Dtheta,z; dXz=-Drz-I1 Q1 Dz(Dzz+dp); dp from all full covariance entries',
        signed='selected principal Bil(W0c,L) is current Xold; remaining targets equal minus selected covariance error minus full axial pressure transport',
        moving='Qe Dz H=Dz(Qe H)+epsilon sigma_e,Z int R^e H; whole partial_Z(c_rho dP) retained',
        fast='zero Haar ACTUAL increments remove only linear Haar terms; dX order min(mu+H+1,2H+1), defects change at min(mu+H,2H)',
        slow='CURRENT slow increments of orders C and C+1 give dp order C but target change C+1, including time viscosity axial pressure and squares',
        target_cycle='C>=C0=2-2kappa maps to C+1/2-kappa; affine slopes and boundary values prove whole half line not stage scan',
        seed='accepted initial target 2-3kappa; real prefix then first sourced C1=3/2-kappa and generic correction give current target C0',
        tails='linear nonzero angular commutators average to zero, complete cut quadratic velocities remain in covariance; base and mean tails not erased or forward inverted',
        domain='unchanged corrected initial tail including ell^6>=8*Cfreeze/gamma0 and frozen-margin Lreq; finite jet and stage constants retained',
        scope='conditional higher-order target bridge; not independent initial field/inverse proof, entire residual audit, old first-band certification or unconditional blow-up')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'current_target_moment_audit';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('cycle_current_target_parent',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():
    m=parent_module();m.validate(json.loads((HERE.parent/'current_target_moment_audit/evidence.json').read_text()))
    for _ in range(12):
        if hasattr(m,'provenance_tools'):break
        m=m.parent_module()
    else:raise ValueError('Missing historical provenance route')
    return m.provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','target_cycle.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'current_target_moment_audit/target_algebra.py',
              HERE.parent/'canonical_initial_mean_join/real_covariance.py',HERE.parent/'complex_signed_covariance_update/covariance.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
@lru_cache(None)
def reviews():return review()
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g=dict(r['checks'])
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Current-cycle record changed')
    if not all(expected['gates'].values()):raise ValueError('Current-cycle gate failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen historical evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates; {len(r['reviews']['current_equation_cancellations']['identities'])} exact current-update identities")
    print('Conditional target recurrence only; retained initial maps and complete residual/forcing audit remain.')
if __name__=='__main__':main()
