#!/usr/bin/env python3
"""Freeze the joint-divisor bound and scoped current centrifugal input split."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_compact_source_quotient','core_compact_input_sharpness','core_fast_inverse_budget',
              'core_next_sourced_mean_tail_join','compact_mean_foundation')
STATUS='mixed-compact-temporal-inverse-and-current-centrifugal-pressure-input-split'
FLAGS={
    'original_joint_integer_divisor_identity_and_same_mode_requirement_proved':True,
    'mixed_raw_source_rule_m_plus_two_p_plus_two_derived_with_full_profiles':True,
    'actual_f0_over_R_all_raw_profile_jets_have_explicit_analytic_recipe':True,
    'five_derivative_saving_only_relative_to_existing_sequential_upper_rules':True,
    'fast_map_compliant_finite_jet_comparison_requires_m_plus_two_p_minus_one':True,
    'full_current_slow_fast_centrifugal_split_and_Haar_variance_retained':True,
    'current_Gtheta_and_rest_programs_derived_from_coarse_inputs_not_fine_callbacks':True,
    'same_current_K0_pressure_moving_rho_and_all_old_domains_retained':True,
    'actual_all_scale_current_source_high_jet_or_support_inputs_validated':False,
    'mixed_finite_jet_comparison_is_canonical_or_uniform_all_jet_counterexample':False,
    'total_original_root_or_normal_DAG_minimal_jet_count_proved':False,
    'whole_A0_pressure_regular_E_or_parent_direct_A1_order_improved':False,
    'next_signed_cycle_absolute_decay_infinite_field_or_blowup_proved':False,
    'historical_evidence_actual_fields_width_h_carriers_cutoffs_or_band_changed':False
}
REMAINING=[
    'validate actual all-band Gtheta and Rrest high-jet weighted/raw norms, current field/row identification and whole-hull support/extension with all fixed profile/base bounds',
    'if using an explicit sum of old prefast rows for Fsum, prove its original cover/label/c normalization transport and norm costs rather than adding differently normalized callbacks',
    'join reduced radial source and complete pressure/support identities before the next generic signed cycle, retaining every older A0/root/normal graph and additive tail',
    'obtain sufficient full regular orders and controlled infinite diagonal, localization/forcing regularity and full velocity growth before any blow-up claim'
]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks
def sources():
    paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_compact_source_quotient/bounds.py',HERE.parent/'core_compact_input_sharpness/bounds.py',
              HERE.parent/'core_fast_inverse_budget/proof_checks.py',HERE.parent/'core_next_sourced_mean_tail_join/budgets.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {k:canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_joint_divisor_profile_nonlinear_pressure_and_budget_controls_pass=bool(r['checks']) and all(r['checks'].values()),
            five_historical_records_bound=len(bound)==5,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Joint divisors, raw profiles, current centrifugal pressure split, scope or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
