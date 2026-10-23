#!/usr/bin/env python3
"""Freeze native transport and the unchanged-band sufficient floor."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_finite_support_closure','core_tangential_normal_closure',
              'core_physical_representative_closure','core_background_finite_head_budget',
              'physical_pulse_foundation','endpoint_foundation')
STATUS='native-axial-transport-retained-power-floor-and-original-slow-gate'
FLAGS={
    'exact_implicit_chart_Euler_axial_derivative_derived':True,
    'actual_heat_exterior_direction_cancels_exactly':True,
    'leading_native_two_profile_first_jet_route_and_same_full_diagonal_retained':True,
    'normal_error_keeps_rounding_and_correction_at_distinct_S_powers':True,
    'original_scalar_band_bound_supplies_uniform_slow_gate':True,
    'actual_profile_suprema_diagonal_identity_or_same_band_floor_certified':False,
    'nearest_nonzero_carrier_replaced_or_new_scale_selected':False,
    'full_pressure_memory_or_radial_cutoff_curl_dropped':False,
    'complete_phase_cone_or_regular_residual_order_promoted':False,
    'physical_absolute_decay_infinite_field_or_blowup_proved':False,
    'historical_frozen_evidence_fields_or_width_changed':False
}
REMAINING=[
    'enclose actual leading native E0/U0 value/X/eta first jets over the original physical hull and retained reference growth/shear/tilt gaps',
    'identify the SAME full low-order Euler diagonal and existing first-band scalar prescription, then evaluate the retained-power normal margin at that supplied band',
    'retain separate same-band target-freezing/root/cone gaps, full radial phase/pressure/cutoff estimates and actual weighted current target/source norms',
    'complete regular physical residual decay, controlled infinite field, localization/forcing regularity and complete velocity growth before blow-up claims'
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
    paths += [HERE.parent/'core_tangential_normal_closure/budgets.py',
              HERE.parent/'core_background_finite_head_budget/budgets.py',
              HERE.parent/'core_physical_representative_closure/budgets.py',
              HERE.parent/'core_finite_support_closure/bounds.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {str(k):canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_chart_transport_and_floor_checks_pass=bool(r['checks']) and all(r['checks'].values()),
            six_historical_records_bound=len(bound)==6,
            all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Native transport, same-band margin, original slow bound, scope or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
