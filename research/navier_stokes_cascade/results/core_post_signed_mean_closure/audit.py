#!/usr/bin/env python3
"""Freeze the conditional post-signed mean block and complete request domains."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_generic_signed_covariance_join','core_preserved_haar_scalar_join','core_final_mean_row_closure',
              'core_fast_row_gain_join','core_fast_inverse_budget','core_coupled_defect_budget',
              'core_reweighted_fresh_mean_closure','core_multiharmonic_mean_residual_budget')
STATUS='conditional-post-signed-fresh-mean-full-regular-row-and-complete-request-domain-closure'
FLAGS={
    'historical_qualitative_finite_cycle_not_presented_as_new':True,
    'new_fast_source_derived_from_old_final_rows_and_complete_signed_covariance':True,
    'actual_A1_velocities_and_positive_exact_time_remainder_retained':True,
    'full_old_fast_correlations_and_three_fresh_current_scalars_kept':True,
    'same_whole_reserved_hull_G_zero_contract_and_preserved_Haar_sum_required':True,
    'all_three_full_regular_mean_rows_conditionally_at_2_minus_2kappa':True,
    'whole_pressure_reconstruction_viscosity_time_and_interactions_kept':True,
    'A1_and_A0_raw_source_requests_and_original_kernel_factors_explicit':True,
    'one_common_pole_covers_entire_new_request_graph_and_older_DAG_references':True,
    'new_2q_4q_operator_constants_are_original_deferred_raw_profile_recipes':True,
    'mean_feedback_order_D_uses_derived_fast_increment_not_coarse_relabelling':True,
    'low_output_signed_pole_used_for_all_huge_compact_source_requests':False,
    'actual_original_weighted_high_jet_profile_hull_or_background_inputs_verified':False,
    'complete_next_nonzero_or_physical_residual_and_infinite_field_certified':False,
    'absolute_decay_unconditional_pulses_or_blowup_proved':False,
    'actual_fields_roots_width_band_or_historical_frozen_evidence_changed':False,
}
REMAINING=[
    'validate actual same-family weighted coefficient/profile/support/background inequalities through the ENTIRE new and older nested request DAGs, with original p chosen for all finite requested jets',
    'join the complete new signed linear pressure/multiharmonic nonzero residual and all physical compact/Gaussian/base tails on the expanded common domains before claiming the full canonical seed',
    'advance common-family regular orders enough for absolute physical decay; the mean block alone does not prove complete residual improvement',
    'join controlled infinite diagonal and localization products, smooth admissible forcing and complete velocity growth before any blow-up claim',
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
    paths += [HERE.parent/'core_fast_row_gain_join/budgets.py',HERE.parent/'core_fast_inverse_budget/proof_checks.py',
              HERE.parent/'core_final_mean_row_closure/budgets.py',HERE.parent/'core_preserved_haar_scalar_join/budgets.py',
              HERE.parent/'core_generic_signed_covariance_join/budgets.py',HERE.parent/'core_actual_scalar_gain_join/budgets.py',
              HERE.parent/'core_reweighted_fresh_mean_closure/budgets.py',HERE.parent/'current_target_moment_audit/target_algebra.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_full_row_scalar_compact_request_and_domain_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   eight_historical_records_bound=len(bound)==8,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Full post-signed mean rows, current scalar/pressure recipe, whole request domains, limits or provenance changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")


if __name__=='__main__':main()
