#!/usr/bin/env python3
"""Freeze the current sourced-block mean/radial and physical-tail join."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_next_sourced_wave_join','core_post_signed_mean_closure','core_post_signed_residual_join',
              'core_preserved_haar_scalar_join','core_final_mean_row_closure','core_fast_inverse_budget',
              'core_physical_tail_join','core_multiharmonic_mean_residual_budget')
STATUS='conditional-next-sourced-block-fresh-means-current-radial-and-physical-tail-join'
FLAGS={
    'actual_new_fast_E_from_derived_postwave_rows_not_fine_input_callback':True,
    'all_three_current_postfast_scalars_include_full_covariance_and_old_fast_correlations':True,
    'Haar_seed_and_all_three_historical_slow_summands_preserved':True,
    'same_current_five_map_whole_G_zero_potential_hull_contract_retained':True,
    'complete_postslow_scalars_C0_plus_one_and_whole_axial_pressure_C0_plus_two':True,
    'all_three_full_regular_mean_rows_conditionally_C0_with_new_A0_source':True,
    'complete_regular_source_E_including_H4_actual_mean_feedback':True,
    'normalized_A1_two_integrations_do_not_erase_old_high_jet_DAGs':True,
    'entire_normalized_and_physical_source_normal_older_target_graphs_traced':True,
    'original_6q_12q_operators_and_all_archived_shifts_and_separate_Ie_retained':True,
    'current_compact_Gaussian_and_assigned_base_tail_allocation_retained':True,
    'same_domain_physical_tail_estimates_require_specific_actual_inputs':True,
    'actual_all_scale_weighted_high_jet_support_background_base_inputs_validated':False,
    'next_generic_signed_cycle_or_infinite_field_supplied':False,
    'complete_actual_residual_improvement_absolute_decay_unconditional_pulses_or_blowup_proved':False,
    'historical_evidence_actual_fields_width_h_roots_or_accepted_band_changed':False,
}
REMAINING=[
    'independently validate actual same-family full-axis/profile/root/frame/background/support/endpoint and assigned base/exterior inequalities at every requested physical scale and nested finite jet',
    'advance the next current generic signed correction and its fresh means on the expanded domains before claiming subsequent regular orders',
    'obtain higher complete regular orders satisfying h*alpha>N+Lambda+2m+eta; present E gives a negative physical exponent, not an actual residual lower bound',
    'join controlled infinite diagonal, global localization products, smooth forcing and complete velocity growth before any blow-up claim',
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
    paths += [HERE.parent/'core_next_sourced_wave_join/budgets.py',HERE.parent/'core_post_signed_mean_closure/budgets.py',
              HERE.parent/'core_post_signed_residual_join/bounds.py',HERE.parent/'core_preserved_haar_scalar_join/budgets.py',
              HERE.parent/'core_preserved_haar_scalar_join/controls.py',HERE.parent/'core_final_mean_row_closure/budgets.py',
              HERE.parent/'core_final_mean_row_closure/controls.py',HERE.parent/'core_fast_row_gain_join/budgets.py',
              HERE.parent/'core_physical_tail_join/bounds.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}

def canonical(value):
    if isinstance(value,Fraction):return str(value)
    if isinstance(value,dict):return {k:canonical(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [canonical(v) for v in value]
    return value

def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_current_mean_scalar_radial_physical_and_request_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   eight_historical_records_bound=len(bound)==8,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Current sourced-block means, pressure/A0, physical graph, scope or provenance changed')

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
