#!/usr/bin/env python3
"""Freeze the conditional finite-support closure without certifying its inputs."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_compact_trace_gate','source_envelope_edge_audit','core_next_sourced_wave_join',
              'compact_mean_foundation','auxiliary_geometry_audit','core_axial_potential_pressure_gain')
STATUS='finite-smooth-support-closure-relative-to-original-weighted-input-contract'
FLAGS={
    'finite_grammar_all_mixed_trace_closure_derived_from_retained_input_contract':True,
    'reciprocal_root_used_only_in_actual_full_weighted_signed_product':True,
    'signed_support_weight_threshold_theta_greater_than_one_half_is_sharp':True,
    'whole_path_source_support_not_replaced_by_entrance_flatness':True,
    'full_cut_potential_curl_and_source_normal_pressure_retained':True,
    'shifted_mean_full_hull_and_global_physical_cover_representatives_retained':True,
    'all_seam_traces_conditional_on_global_input_contract_not_independent_samples':True,
    'original_weighted_canonical_input_or_same_band_floor_gaps_certified':False,
    'old_raw_norm_request_graph_or_regular_source_order_improved':False,
    'raw_outgoing_forced_response_assumed_compact':False,
    'negative_physical_time_zero_extension_assumed':False,
    'absolute_residual_decay_infinite_field_or_blowup_proved':False,
    'historical_evidence_actual_fields_or_first_band_changed':False
}
REMAINING=[
    'validate weighted mixed target/source inequalities at the earliest unvalidated actual node of the retained finite graph, with original canonical root/frame data',
    'check original unchanged-band normal/root/cone/freezing and whole closed-path/copy compatibility gates; no later first band substitutes for them',
    'retain all huge raw norm/operator/profile requests and complete pressure/mean/Haar histories; smooth trace closure alone does not bound their constants',
    'then advance the signed cycle, complete physical residual decay, infinite field, localization/force regularity and complete velocity growth before blow-up claims'
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
    paths += [HERE.parent/'source_envelope_edge_audit/source_bounds.py',
              HERE.parent/'core_next_sourced_wave_join/budgets.py',
              HERE.parent/'core_axial_potential_pressure_gain/bounds.py',
              HERE.parent/'core_compact_trace_gate/bounds.py',
              HERE.parent/'auxiliary_geometry_audit/geometry.py',
              HERE.parent/'compact_mean_corrections/mean_system.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {str(k):canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_factor_closure_and_scope_checks_pass=bool(r['checks']) and all(r['checks'].values()),
            six_historical_records_bound=len(bound)==6,
            all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Finite support implication, original input gate, pressure or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
