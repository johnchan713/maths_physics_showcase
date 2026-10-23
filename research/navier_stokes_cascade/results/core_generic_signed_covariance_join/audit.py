#!/usr/bin/env python3
"""Freeze current generic signed covariance, pressure domains and target gain."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_preserved_haar_scalar_join','core_final_mean_row_closure','core_primary_root_budget',
              'core_signed_target_edge_budget','core_reweighted_fresh_mean_closure','core_localized_stress_assembly',
              'current_target_cycle_bridge','core_multiharmonic_mean_residual_budget')
STATUS='conditional-current-generic-signed-covariance-pressure-domain-and-target-join'
FLAGS={
    'historical_qualitative_generic_target_gain_not_presented_as_new':True,
    'original_root_localized_principal_cross_identity_with_cutoffs_and_partition_kept':True,
    'all_six_covariance_entries_complex_crosses_and_signed_squares_retained':True,
    'same_current_final_row_and_whole_moment_recipes_supply_target_inputs':True,
    'both_distinct_target_components_and_their_edge_allowances_reach_original_root':True,
    'conditional_post_signed_target_and_Haar_tangential_order_is_2_minus_2kappa':True,
    'same_moving_full_pressure_and_projection_derivatives_retained':True,
    'new_pressure_and_target_primitive_domains_supplied_by_original_kernel_recipes':True,
    'scalar_defect_radial_Haar_and_full_pointwise_pressure_not_promoted':True,
    'original_primary_or_error_bound_inferred_from_updated_wave_norm':False,
    'new_actual_fields_roots_width_supports_or_accepted_band_introduced':False,
    'actual_weighted_target_root_common_chart_high_jet_and_support_inputs_verified':False,
    'subsequent_fresh_fast_slow_means_or_complete_cycle_certified':False,
    'absolute_residual_decay_unconditional_pulses_infinite_field_or_blowup_proved':False,
    'historical_frozen_evidence_changed':False,
}
REMAINING=[
    'validate same-family actual weighted current target, original principal/error summands, complete signed/root/common-chart coefficient and compact endpoint/support inputs at all requested jets',
    'derive actual post-signed fast then fresh slow difference bounds and all three scalar updates, radial A0 repair and new compact-tail source domains before claiming the complete C0 canonical seed',
    'establish sufficiently high common-family regular orders and complete physical residual estimates; current source order D does not certify absolute physical decay',
    'join controlled infinite diagonal and localization products, forcing regularity and complete velocity growth before any blow-up claim',
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
    paths += [HERE.parent/'core_primary_root_budget/budgets.py',HERE.parent/'core_signed_target_edge_budget/budgets.py',
              HERE.parent/'core_reweighted_fresh_mean_closure/budgets.py',HERE.parent/'core_preserved_haar_scalar_join/budgets.py',
              HERE.parent/'core_final_mean_row_closure/budgets.py',HERE.parent/'current_target_cycle_bridge/target_cycle.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_signed_covariance_pressure_domain_and_current_target_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   eight_historical_records_bound=len(bound)==8,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Generic signed covariance, pressure domains, target gain, limits or provenance changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")


if __name__=='__main__':main()
