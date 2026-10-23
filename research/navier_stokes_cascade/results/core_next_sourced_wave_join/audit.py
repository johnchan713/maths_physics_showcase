#!/usr/bin/env python3
"""Freeze the next current-source constrained wave and expanded domains."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_post_signed_residual_join','core_post_signed_mean_closure','core_first_sourced_inverse_budget',
              'core_tangent_propagator_closure','core_sourced_normal_transfer','core_sourced_covariance_moment_join',
              'core_multiharmonic_mean_residual_budget','sourced_inverse_foundation')
STATUS='conditional-next-current-source-eight-mode-constrained-wave-and-wave-stage-gain'
FLAGS={
    'same_derived_current_D_source_program_not_new_fine_source_assumption':True,
    'eight_signed_input_harmonics_original_inverse_and_full_source_normal_pressure':True,
    'Gaussian_adapted_ambient_jets_use_full_length_not_high_mode_residence':True,
    'same_cut_potential_and_pressure_full_graph_curl_and_all_six_covariances':True,
    'complete_velocity_D_distinct_from_pressure_potential_C0':True,
    'wave_stage_regular_gain_D_to_2_minus_3kappa_before_fresh_means':True,
    'full_tangential_E_Haar_target_scalar_C0_and_whole_axial_pressure_kept':True,
    'H4_ratio4_and_eighth_output_not_old_single_generated_mode_list':True,
    'new_6q_pressure_domain_and_6q_12q_kernel_recipes_plus_all_old_shifts':True,
    'whole_path_and_flat_entrance_support_contracts_explicit':True,
    'new_A0_pressure_residual_not_replaced_by_stale_prewave_tail':True,
    'actual_current_source_high_jet_or_all_scale_support_inputs_instantiated':False,
    'all_three_new_full_pointwise_mean_rows_or_radial_A0_bound_joined':False,
    'subsequent_fresh_means_or_generic_signed_full_cycle_complete':False,
    'complete_actual_residual_improvement_absolute_decay_infinite_field_or_blowup_proved':False,
    'actual_band_roots_geometry_width_or_historical_evidence_changed':False,
}
REMAINING=[
    'validate actual same-family full-path background/frame/root/coefficient/profile and all nested high-jet/support/flat-entrance inputs; current norm DAG remains deferred',
    'apply actual fast map to derived postwave Qcirc tangential E rows, recompute current postfast pressure and all three scalars, then apply reserved five-map with whole-hull G=0 contract',
    'join NEW pointwise radial A0 and physical compact/Gaussian/base tails with full expanded request graphs before claiming the next whole residual/cycle',
    'derive later generic signed cycle and common-family higher regular orders, controlled infinite diagonal, localization/forcing regularity and complete velocity growth before any blow-up claim',
]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks

def sources():
    paths=[HERE/n for n in ('README.md','budgets.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_post_signed_residual_join/bounds.py',HERE.parent/'core_post_signed_mean_closure/budgets.py',
              HERE.parent/'core_first_sourced_inverse_budget/budgets.py',HERE.parent/'core_tangent_propagator_closure/budgets.py',
              HERE.parent/'core_sourced_normal_transfer/budgets.py',HERE.parent/'core_sourced_covariance_moment_join/budgets.py',
              HERE.parent/'core_multiharmonic_mean_residual_budget/budgets.py',HERE.parent/'finite_induction_foundation/finite_checks.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}

def canonical(value):
    if isinstance(value,Fraction):return str(value)
    if isinstance(value,dict):return {k:canonical(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [canonical(v) for v in value]
    return value

def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_current_inverse_harmonic_covariance_and_domain_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   eight_historical_records_bound=len(bound)==8,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Current sourced inverse, full pressure/curl, next domains, scope or provenance changed')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    payload=json.dumps(r,indent=2,sort_keys=True)+'\n'
    with a.output.open('x') as f:f.write(payload)
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")

if __name__=='__main__':main()
