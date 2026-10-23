#!/usr/bin/env python3
"""Freeze finite-prefix physical tails and the unclosed absolute regular gap."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_final_mean_row_closure','core_fast_inverse_budget','current_residual_ledger_audit',
              'physical_pulse_transfer','auxiliary_geometry_audit','core_multiharmonic_mean_residual_budget','core_fast_row_gain_join')
STATUS='conditional-current-prefix-physical-tail-transfer-and-absolute-regular-order-gap'
FLAGS={
    'exact_compact_residual_signs_and_both_angular_projections_retained':True,
    'actual_A1_velocities_stay_in_regular_nonlinear_interactions':True,
    'original_temporal_Gaussian_ledger_is_additive_and_never_inverted':True,
    'all_fixed_physical_compact_tail_recipes_include_normalization_and_graph_losses':True,
    'same_domain_polynomial_and_band_to_natural_log_Gaussian_constants_derived':True,
    'source_raw_profiles_original_poles_width_gap_and_nested_high_jets_retained':True,
    'physical_regular_upper_bound_has_explicit_negative_current_exponent':True,
    'precise_sufficient_absolute_regular_order_inequality_identified':True,
    'spatial_diagonal_and_global_localization_cutoffs_keep_separate_estimates':True,
    'new_analytical_assumption_or_derivative_dependent_band_added':False,
    'illustrative_h_or_huge_jet_counts_are_actual_input_evaluations':False,
    'actual_residual_lower_bound_or_route_refutation_proved':False,
    'future_full_induction_arbitrary_diagonal_or_infinite_field_supplied':False,
    'actual_profile_support_base_exterior_and_full_physical_inputs_verified':False,
    'absolute_residual_decay_complete_improvement_unconditional_pulses_or_blowup_proved':False,
    'historical_frozen_evidence_fields_or_width_changed':False,
}
REMAINING=[
    'independently validate actual background/cone/root/profile, all requested high-jet and support/endpoint inputs and the exact assigned base/exterior flat remainder',
    'establish common-family later regular orders meeting h*alpha > N+Lambda+2m+eta; current D does not certify physical decay and the exponent bound is not a residual lower bound',
    'join controlled shrinking-diagonal and force-localization cutoff product jets to one infinite field; stage-dependent constants do not permit an arbitrary diagonal',
    'prove forcing regularity and complete velocity growth before any blow-up claim',
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
    paths += [HERE.parent/'core_final_mean_row_closure/budgets.py',HERE.parent/'core_fast_inverse_budget/proof_checks.py',
              HERE.parent/'current_residual_ledger_audit/residual_ledger.py',HERE.parent/'core_fast_row_gain_join/budgets.py',
              HERE.parent/'physical_pulse_foundation/README.md',HERE.parent/'geometric_pulse_replay/README.md']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_physical_tail_allocation_exponent_and_constant_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   seven_historical_records_bound=len(bound)==7,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Physical tail allocation, constants, regular-order gap, input limits, provenance or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")


if __name__=='__main__':main()
