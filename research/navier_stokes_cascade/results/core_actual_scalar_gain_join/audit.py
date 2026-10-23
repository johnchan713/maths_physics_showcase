#!/usr/bin/env python3
"""Freeze the derived actual scalar gain and its retained support limitation."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_multiharmonic_mean_residual_budget','core_reweighted_fresh_mean_closure',
              'core_sourced_covariance_moment_join','core_signed_prefix_source_budget',
              'core_coupled_defect_budget','core_dynamic_mean_recheck','compact_mean_foundation')
STATUS='derived-current-scalar-C-from-existing-coarse-actual-means-with-retained-reserved-hull-contract'
FLAGS={
    'all_three_postslow_scalar_remainders_at_two_B_derived':True,
    'subunit_B_plus_one_promotion_rejected':True,
    'postfast_scalar_C_composed_without_new_incoming_scalar_callback':True,
    'fresh_slow_C_size_follows_from_same_actual_five_map':True,
    'complete_axial_velocity_covariance_and_weighted_pressure_moment_retained':True,
    'whole_moving_compensator_product_and_actual_fast_A1_fields_retained':True,
    'entire_reserved_potential_hull_G_zero_is_retained_original_premise':True,
    'nonzero_G_support_violation_has_linear_scalar_leaks':True,
    'higher_pole_moment_recipes_and_long_shell_embeddings_retained':True,
    'new_analytical_assumption_or_accepted_band_gate_added':False,
    'original_actual_support_or_input_enclosures_independently_certified':False,
    'actual_fast_row_D_gain_or_combined_feedback_D_proved':False,
    'pointwise_mean_pressure_C_or_complete_residual_gain_proved':False,
    'unconditional_pulses_infinite_field_or_blowup_proved':False,
    'historical_evidence_fields_width_or_supports_changed':False,
}
REMAINING=[
    'prove actual post-wave zero-Haar tangential mean-row inputs at D or the parent exact combined feedback inequality (8); fresh fast feedback still has only B',
    'validate original same-field primitive root/frame coefficient and recipe inputs and the entire reserved potential-hull support/G-zero/power-law contracts on every physical scale',
    'join all pointwise mean equations and full pressure to the complete residual including viscosity, spatial cutoffs, exact additive tails and every interaction without forwarding tails to an inverse',
    'establish uniform actual pulse family, infinite construction and forcing/blow-up obligations only after their complete estimates',
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
    paths += [HERE.parent/n/'proof_checks.py' for n in ('core_coupled_defect_budget','core_dynamic_mean_recheck','core_slow_inverse_budget')]
    paths += [HERE.parent/'core_reweighted_fresh_mean_closure/controls.py',HERE.parent/'complete_covariance_initial_means/means.py',
              HERE.parent/'compact_mean_foundation/mean_maps.py',HERE.parent/'complex_mean_relaxation/README.md']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_subunit_scalar_and_same_recipe_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   seven_historical_records_bound=len(bound)==7,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Actual scalar gain, support limitation, provenance or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")


if __name__=='__main__':main()
