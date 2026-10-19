#!/usr/bin/env python3
"""Explicit repair jets and the precise retained cutoff-shape input gap."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from derivatives import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='exact-repair-bump-jets-and-conditional-cutoff-shape-budget'
DEPENDENCIES=('global_axial_angular_bound','first_seminorm_bound','coefficient_support_foundation',
    'annular_first_correction','matching_root_foundation')
FLAGS=dict(actual_unit_dy_mass_squared_flat_bumps_bounded_through_two_log_derivatives=True,
    disjoint_five_bump_supports_and_reserved_mean_interval_unchanged=True,
    retained_inner_cutoff_shape_derivative_budget_is_a_separate_missing_input=True,
    monotone_compressed_family_disproves_uniform_upper_bound_from_support_alone=True,
    cutoff_log_Jacobian_and_E_to_F_inverse_radius_terms_retained=True,
    same_first_jet_interface_composed_conditionally_on_actual_shape_norms=True,
    complete_pressure_source_rows_gauge_and_both_mean_memories_preserved=True,
    no_actual_field_root_frequency_width_cutoff_or_support_replaced=True,
    historical_frozen_evidence_preserved=True,
    canonical_sigma_shape_identified_with_unchanged_inner_cutoff=False,
    missing_shape_budget_is_nonexistence_of_suitable_cutoff=False,
    actual_cutoff_derivative_enclosures_or_combined_scale_gate_verified=False,
    actual_giant_parameter_profile_numerically_materialized=False,
    complete_first_stress_leading_pointwise_and_shifted_axial_inputs_closed=False,
    unconditional_pulses_complete_residual_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['identify_or_bound_actual_retained_cutoff_normalized_shape_L1_L2_without_silently_replacing_it',
    'bound_remaining_actual_pointwise_w0_shifted_axial_and_leading_radial_inputs_then_join_first_Cstar',
    'verify_unchanged_tilt_and_combined_scale_freezing_geometry_inverse_and_acceptance_gates',
    'discharge_uniform_actual_pulse_complete_covariance_mean_and_viscous_residual_obligations']
def protocol():
    return dict(parent_commit='ec57a0079d884edb2b1be30a587f98455628b1ef',
        cutoff='same Xkeep=Xa exp(tc/8),Xcut=Xa exp(3tc/16),tc=C^-20; normalized log width Delta=tc/16',
        shape='L1=sup|K_prime|,L2=sup|K_second|; necessary L1>=1,L2>=4 are not upper bounds',
        obstruction='K_n=1-sigma(n(s-1/2)+1/2), norms >=8n and >=16n^2; auxiliary family, no actual cutoff replaced',
        actual_step='squared-flat logistic, derivatives <=(32,1216,220000) proved by global real-variable envelopes',
        actual_bumps='same five dy-unit-mass width 1/10 bumps; (B0,B1,B2)<=(320,121600,220000000)',
        composition='same inverse 800 and first seminorm extended_jets; rho3 lower radius, r=sqrt(rho3); pressure reconstructed, never cut',
        scope='conditional first-jet input and precise cutoff-shape gap; no same-field scale, full residual or blow-up certificate')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('cutoff_repair_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','derivatives.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py',HERE.parent/'first_seminorm_bound/seminorms.py',
        HERE.parent/'stress_realization_audit/repair.py',HERE.parent/'outer_pressure_pilot/schedule.py',
        HERE.parent/'annular_first_correction/extension.py']
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Cutoff/repair derivative record changed')
    if not all(expected['gates'].values()):raise ValueError('Cutoff/repair derivative check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual cutoff shape and combined first-stress scale remain unverified.')
if __name__=='__main__':main()
