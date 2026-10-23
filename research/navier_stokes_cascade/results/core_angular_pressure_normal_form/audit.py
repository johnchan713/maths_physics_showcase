#!/usr/bin/env python3
"""Freeze the current angular normal form and scoped full mean-pressure bound."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_mixed_compact_time_budget','core_compact_source_quotient',
              'core_next_sourced_mean_tail_join','core_preserved_haar_scalar_join','core_fast_inverse_budget')
STATUS='current-angular-row-normal-form-and-conditional-order-one-mean-pressure'
FLAGS={
    'final_angular_row_yields_exact_full_accumulated_swirl_normal_form':True,
    'variable_radial_multiplier_derivative_and_all_graph_viscosities_retained':True,
    'new_local_pressure_K1_is_restored_not_annihilated_with_A0':True,
    'Haar_variance_moving_moment_and_complete_old_pressure_memory_retained':True,
    'order_one_current_mean_pressure_proved_conditionally_on_existing_inputs':True,
    'fine_angular_input_is_existing_pressure_free_row_not_current_radial_output':True,
    'four_distinct_compact_source_requests_and_all_old_domains_retained':True,
    'actual_all_band_high_jet_Haar_profile_or_support_norms_validated':False,
    'total_original_minimal_jet_or_full_residual_order_improved':False,
    'total_pressure_or_parent_direct_A1_order_promoted':False,
    'absolute_decay_next_signed_cycle_infinite_field_or_blowup_proved':False,
    'historical_frozen_evidence_actual_fields_width_band_or_roots_changed':False
}
REMAINING=[
    'validate same-field actual current J,F,v^2/R,Rrest bounds at the distinct raw orders, including derived angular-row and sharper preserved Haar inputs',
    'validate actual V_R kernel jets and every whole compact support/trace, original weighted temporal and pressure/profile conversion on the unchanged accepted band',
    'retain complete K1 local pressure, moving rho P, covariance and axial pressure memory when joining the next generic signed cycle',
    'obtain sufficient complete regular orders, controlled infinite diagonal, localization/forcing regularity and full velocity growth before absolute decay or blow-up claims'
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
    paths += [HERE.parent/'core_mixed_compact_time_budget/bounds.py',HERE.parent/'core_next_sourced_mean_tail_join/budgets.py',
              HERE.parent/'core_preserved_haar_scalar_join/budgets.py',HERE.parent/'compact_mean_corrections/mean_system.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {k:canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_full_row_pressure_profile_and_budget_controls_pass=bool(r['checks']) and all(r['checks'].values()),
            five_historical_records_bound=len(bound)==5,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Angular normal form, full local pressure, input scope or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
