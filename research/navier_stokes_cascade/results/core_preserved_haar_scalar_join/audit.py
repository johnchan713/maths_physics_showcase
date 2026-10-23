#!/usr/bin/env python3
"""Freeze the complete preserved-Haar scalar join and its actual input gap."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_physical_tail_join','core_final_mean_row_closure','core_actual_scalar_gain_join',
              'core_initial_scalar_seed_budget','core_scalar_cycle_budget','core_coupled_defect_budget',
              'core_haar_row_stress_join')
STATUS='conditional-complete-preserved-Haar-current-scalar-refinement-and-weighted-input-gap'
FLAGS={
    'known_qualitative_Haar_preservation_not_presented_as_new':True,
    'all_three_complete_current_postslow_scalars_have_conditional_C_plus_one_budget':True,
    'old_Haar_norm_is_derived_finite_slow_sum_not_coarse_full_norm_relabelled':True,
    'first_signed_complete_covariance_and_full_old_fast_correlations_retained':True,
    'current_fresh_five_map_and_latest_scalar_derived_slow_increment_used':True,
    'axial_velocity_square_and_whole_moving_pressure_moment_kept':True,
    'original_seed_signed_weight_transfer_and_support_input_inequalities_explicit':True,
    'original_high_jet_compact_tail_domains_and_regular_D_unchanged':True,
    'new_fine_output_bound_assumed_or_new_band_imposed':False,
    'original_weighted_seed_signed_inputs_or_hull_independently_verified':False,
    'total_pressure_or_fast_velocity_promoted_to_sharp_Haar_order':False,
    'next_complete_signed_cycle_or_full_physical_residual_closed':False,
    'unconditional_pulse_absolute_decay_infinite_field_or_blowup_proved':False,
    'historical_frozen_evidence_or_actual_fields_changed':False,
}
REMAINING=[
    'validate actual same-family weighted initial SeedBudget fields/defects, complete cut original-root signed coefficients and all chart/weight/hull transfer inequalities at every required finite jet',
    'compose the generic original-root signed update with current complete targets and whole moments; prove complete covariance cancellation at 2-2kappa, then fresh mean rows and compact-tail bounds',
    'close common-family later regular orders and actual high-jet/support/base inputs needed for absolute physical decay; current nonzero regular order remains D',
    'join controlled physical diagonal and localization cutoff products, forcing regularity and velocity growth before any blow-up claim',
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
    paths += [HERE.parent/'core_actual_scalar_gain_join/budgets.py',HERE.parent/'core_initial_scalar_seed_budget/proof_checks.py',
              HERE.parent/'core_final_mean_row_closure/budgets.py',HERE.parent/'core_coupled_defect_budget/proof_checks.py',
              HERE.parent/'core_signed_target_edge_budget/README.md',HERE.parent/'core_signed_prefix_source_budget/budgets.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_complete_scalar_Haar_pressure_rate_and_composition_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   seven_historical_records_bound=len(bound)==7,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Complete scalar identity, finite Haar recipe, weighted input gap, scope or provenance changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")


if __name__=='__main__':main()
