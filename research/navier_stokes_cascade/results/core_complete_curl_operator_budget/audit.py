#!/usr/bin/env python3
"""Freeze the complete graph/cylindrical curl and consistent radial-chart budgets."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_canonical_haar_column_jet_closure','core_canonical_entrance_jet_closure',
    'core_reference_scaled_normal_floor','core_phase_coefficient_budget','core_primary_root_budget',
    'auxiliary_geometry_audit','canonical_initial_mean_join','core_primary_aggregate_atlas')
STATUS='conditional-complete-curl-operator-and-radial-chart-budget'
FLAGS=dict(full_cylindrical_connection_graph_axial_and_potential_derivative_terms_retained=True,
    curl_operator_input_derived_from_full_normal_reciprocal_and_original_scale_identities=True,
    all_coefficient_derivatives_and_extra_normal_jet_in_D0_retained=True,
    U_chart_conjugates_first_order_operator_before_composing_every_root_input=True,
    same_physical_fields_roots_pressure_memory_and_supports_retained=True,
    full_background_n_plus_thirteen_whole_axial_moment_ledger_explicit=True,
    all_six_covariances_axial_mean_and_reserved_pressure_support_obligations_retained=True,
    actual_background_suprema_weighted_target_or_cutoff_input_certified=False,
    original_accepted_band_length_carriers_tilt_width_or_historical_evidence_changed=False,
    unconditional_pulse_complete_residual_infinite_field_or_blowup_proved=False)
REMAINING=[
    'derive actual weighted target/loghalf and fixed cutoff jet inputs with original value/cone hypotheses, or identify their precise blocking inequality',
    'retain validity of full background/phase primitives and original same physical accepted domains, with the explicit higher derivative ledger',
    'join complete-wave covariance and fresh means including axial velocity, Pi1/Pi2 pressure memory, moving moments and reserved supports',
    'close every viscosity/cutoff/interaction residual and infinite realization before unconditional pulse bounds, full physical residual improvement or blow-up claims']
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
    paths += [HERE.parent/p for p in ('core_canonical_haar_column_jet_closure/bounds.py',
        'core_canonical_entrance_jet_closure/bounds.py','core_reference_scaled_normal_floor/bounds.py',
        'core_phase_coefficient_budget/budgets.py','core_first_sourced_inverse_budget/budgets.py',
        'core_primary_root_budget/budgets.py','auxiliary_geometry_audit/geometry.py',
        'actual_band_pulse_bounds/README.md')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {str(k):canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_curl_chart_and_join_controls=bool(r['checks']) and all(r['checks'].values()),
            eight_historical_records_bound=len(bound)==8,
            all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Complete cylindrical/graph curl, chart conjugation, root join, ledger or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
