#!/usr/bin/env python3
"""Freeze the complete Haar stress criterion and mixed postslow radial gain."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_fast_row_gain_join','core_actual_scalar_gain_join','core_initial_haar_row_budget',
              'core_sourced_covariance_moment_join','canonical_initial_mean_join','core_primary_root_budget',
              'core_tangent_propagator_closure','compact_mean_foundation')
STATUS='complete-Haar-row-stress-criterion-and-conditional-postslow-radial-Haar-B-plus-C-gain'
FLAGS={
    'complete_Haar_row_envelope_with_all_mean_correlations_and_pressure_derived':True,
    'nongap_Haar_terms_have_order_B_plus_one_without_graph_loss':True,
    'principal_stress_gap_and_fine_flux_C_separated_explicitly':True,
    'same_field_canonical_identity_suffices_conditionally_for_Haar_C':True,
    'two_different_primary_families_need_not_be_identified_for_canonical_route':True,
    'positive_transverse_equal_energy_equal_scalar_countercontrol_has_nonzero_Haar_rows':True,
    'scalar_and_closeness_bounds_alone_do_not_certify_selected_spatial_stress':True,
    'postslow_three_scalars_and_radial_Haar_compensator_B_plus_C_derived':True,
    'whole_moving_pressure_product_and_pointwise_A0_remainder_retained':True,
    'proper_poles_profiles_and_exact_tail_ledger_unchanged':True,
    'new_analytical_assumption_or_accepted_band_gate_added':False,
    'actual_whole_support_principal_stress_identity_or_frame_transfer_certified':False,
    'tensor_countercontrol_is_the_prescribed_canonical_pulse_or_cone':False,
    'pointwise_pressure_or_pointwise_A0_bound_promoted_to_B_plus_C':False,
    'original_actual_input_support_high_jets_or_full_residual_gain_certified':False,
    'unconditional_pulses_infinite_field_or_blowup_proved':False,
    'historical_evidence_fields_width_or_supports_changed':False,
}
REMAINING=[
    'verify principal stress matching and all localization/partition normalizations for the SAME family used in every ledger, or consistently carry out the full canonical rebuilt route; do not transfer another family by closeness',
    'validate original actual all-fixed-order coefficient/profile and support/endpoint contracts including the parents large finite A1 jet dependencies',
    'join complete pointwise A0/A1 and all exact Gaussian/spatial/diagonal cutoff tails to the physical residual without forwarding additive tails to a sourced inverse',
    'prove forcing regularity, uniform actual pulse family, infinite realization and complete-velocity growth before any blow-up claim',
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
    paths += [HERE.parent/n/'budgets.py' for n in ('core_fast_row_gain_join','core_actual_scalar_gain_join',
              'core_initial_haar_row_budget','core_sourced_covariance_moment_join','core_primary_root_budget')]
    paths += [HERE.parent/'core_initial_haar_row_budget/controls.py',HERE.parent/'core_dynamic_mean_recheck/proof_checks.py',
              HERE.parent/'core_reweighted_fresh_mean_closure/controls.py',HERE.parent/'core_reweighted_fresh_mean_closure/budgets.py',
              HERE.parent/'complete_covariance_initial_means/means.py',HERE.parent/'compact_mean_foundation/mean_maps.py',
              HERE.parent/'current_target_cycle_bridge/README.md']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_spatial_stress_mixed_rate_and_composition_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   eight_historical_records_bound=len(bound)==8,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Haar stress criterion, radial gain, actual identity limitation, provenance or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")


if __name__=='__main__':main()
