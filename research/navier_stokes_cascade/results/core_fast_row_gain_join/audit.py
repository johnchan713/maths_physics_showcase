#!/usr/bin/env python3
"""Freeze the conditional actual fast-row and finite regular-source gain join."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_actual_scalar_gain_join','core_multiharmonic_mean_residual_budget',
              'core_reweighted_fresh_mean_closure','core_fast_inverse_budget',
              'core_signed_prefix_source_budget','core_sourced_covariance_moment_join',
              'current_target_cycle_bridge','compact_mean_foundation')
STATUS='conditional-actual-fast-row-D-and-finite-regular-source-gain-with-explicit-A1-high-jets'
FLAGS={
    'full_three_component_pressure_difference_before_estimating_derived':True,
    'actual_fast_time_cancellation_keeps_nonzero_A1_current_source_remainder':True,
    'old_fast_and_slow_differences_at_H_from_coarse_B_fields_derived':True,
    'postwave_zero_Haar_fast_rows_D_and_fresh_actual_fast_D_derived_conditionally':True,
    'parent_scalar_C_and_actual_slow_C_join_retained':True,
    'conditional_mean_feedback_D_and_next_regular_source_D_follow':True,
    'large_finite_A1_source_and_raw_kernel_dependencies_explicit':True,
    'all_pressure_viscosity_cylindrical_connections_and_mean_crosses_retained':True,
    'complete_sourced_covariance_and_proper_operator_poles_retained':True,
    'original_exact_additive_tails_excluded_from_sourced_inverse':True,
    'new_analytical_assumption_or_accepted_band_gate_added':False,
    'high_jet_actual_input_or_raw_kernel_polynomials_numerically_enclosed':False,
    'original_actual_support_or_primary_frame_contract_independently_certified':False,
    'all_zero_angular_equations_and_complete_residual_gain_certified':False,
    'unconditional_pulses_infinite_field_or_blowup_proved':False,
    'historical_evidence_fields_width_or_supports_changed':False,
}
REMAINING=[
    'independently validate original actual coefficient all-fixed-order profile and same-primary/frame contracts, including the large finite jets exposed by the A1 recipe and full support/endpoint extensions',
    'join the full angular-zero Haar target and pressure equations with the conditional nonzero regular gain; scalar moment and Q1 gains alone are not the full mean residual',
    'bound the entire unchanged exact additive tail ledger, all physical derivatives and spatial/diagonal cutoffs and interactions without forwarding additive tails to an inverse',
    'prove uniform actual pulse family, forcing regularity, infinite field and complete-velocity growth before any blow-up claim',
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
    paths += [HERE.parent/n/'budgets.py' for n in ('core_reweighted_fresh_mean_closure','core_sourced_covariance_moment_join','core_signed_prefix_source_budget')]
    paths += [HERE.parent/n/'proof_checks.py' for n in ('core_fast_inverse_budget','core_fixed_profile_tame')]
    paths += [HERE.parent/'core_reweighted_fresh_mean_closure/controls.py',HERE.parent/'complete_covariance_initial_means/means.py',
              HERE.parent/'compact_mean_foundation/mean_maps.py',HERE.parent/'current_target_cycle_bridge/README.md']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_subtraction_A1_rate_and_composition_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   eight_historical_records_bound=len(bound)==8,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Actual fast-row gain, finite A1 dependencies, provenance or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")


if __name__=='__main__':main()
