#!/usr/bin/env python3
"""Current full-covariance target identities and the higher-order cancellation gate."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from target_algebra import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='current-target-full-weight-identities-with-higher-order-cancellation-gate'
DEPENDENCIES=('source_envelope_edge_audit','compact_mean_foundation','complex_mean_relaxation',
              'canonical_initial_mean_join','canonical_finite_cycle_join')
FLAGS=dict(current_full_covariance_moment_map_includes_radial_time_and_axial_transport=True,
    both_current_velocity_moment_constraints_and_whole_moving_pressure_product_retained=True,
    signed_target_exact_radial_flux_plus_compensated_transport_primitive=True,
    all_fixed_jet_full_zeta_functional_membership_relative_to_current_field_bounds=True,
    sharper_epsilon_order_requires_current_residual_cancellation_inequality=True,
    smooth_raw_state_control_blocks_order_gain_from_support_and_moments_alone=True,
    old_new_mean_correlations_actual_A0_A1_and_same_base_Pi1_Pi2_retained=True,
    historical_sources_records_and_corrected_freezing_contract_preserved=True,
    manufactured_symbolic_controls_separate_from_continuum_derivations=True,
    higher_order_current_target_cancellation_established_for_every_canonical_cycle=False,
    raw_state_control_is_counterexample_to_actual_canonical_construction=False,
    historical_first_band_or_force_recertified=False,
    actual_giant_fields_constants_or_first_band_evaluated=False,
    unconditional_pulse_complete_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['prove_current_weighted_Ae_plus_Ie_he_cancellation_at_each_actual_canonical_signed_substep_not_just_raw_field_membership',
    'join_that_estimate_to_complete_current_residual_and_exact_tail_reclassification_with_all_covariances_pressure_viscosity_and_mean_correlations',
    'discharge_retained_analytical_construction_hypotheses_before_any_unconditional_claim']

def protocol():
    return dict(parent_commit='08caced6f10105d98faf59af4e251af26c3a2791',
      current='Haar products averaged after multiplication, six full covariance entries, actual means and same complete base',
      moments='P=tau int B+int(U-F)/R-epsilon partial_Z int H; Qg=tau int R^2 B+int R(F+U)-epsilon partial_Z int R^2 H',
      pressure='compact pbar primitive has zero total by construction; int R pbar=-Qg/2+c_rho P; actual unaveraged A0 retained',
      signed='Xi_e=-Ae-Ie he, with he=epsilon partial_Z He-tau ue-sigma_e Me and exact zero weighted integral',
      inequality='[Xi_e]_m<=C_m([Ae]_m+[he]_m); higher order needs [Ae+Ie he]_m<=C_jm epsilon^alpha_j',
      obstruction='compact smooth swirl control at a(T0)=0,a_prime(T0)=1 has Xi_theta=epsilon^(mu+1) psi/R^2, no higher exponent follows from raw state',
      domain='unchanged corrected initial selector includes ell^6>=8*Cfreeze/gamma0 and frozen margin Lreq; no field or force identification',
      scope='functional target membership and exact current identities, not all later epsilon-order cancellations or unconditional proof')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'source_envelope_edge_audit';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('current_target_edge_parent',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old

@lru_cache(None)
def history():
    m=parent_module();m.validate(json.loads((HERE.parent/'source_envelope_edge_audit/evidence.json').read_text()))
    for _ in range(10):
        if hasattr(m,'provenance_tools'):break
        m=m.parent_module()
    else:raise ValueError('Historical provenance route absent')
    return m.provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','target_algebra.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    paths += [HERE.parent/'compact_mean_corrections/README.md']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Current-target record changed')
    if not all(expected['gates'].values()):raise ValueError('Current-target gate failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; historical records are not overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates; {len(r['reviews']['exact_current_state']['identities'])} exact current identities")
    print('Higher-order canonical target cancellations and unconditional proof remain unverified.')
if __name__=='__main__':main()
