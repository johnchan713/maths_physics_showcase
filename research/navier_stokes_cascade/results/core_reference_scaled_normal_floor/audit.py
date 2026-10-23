#!/usr/bin/env python3
"""Freeze the corrected normal-floor frontier and complete inverse transfer."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_leading_direction_enclosure','core_normal_length_closure',
    'core_sourced_normal_transfer','core_tangent_propagator_closure',
    'core_primary_tail_recheck','core_first_sourced_inverse_budget',
    'frame_column_foundation','core_native_axial_transport_floor')
STATUS='reference-scaled-normal-floor-on-original-frame-tail'
FLAGS=dict(full_original_tangential_error_retains_reference_relative_information=True,
    original_Riccati_slack_gives_31_over_32_reference_floor=True,
    original_carrier_calibration_recovers_phase_beta_bound_without_extra_margin=True,
    stronger_bound_uses_full_error_proof_not_coarse_uniform_floor_alone=True,
    recent_native_margin_is_optional_for_original_full_frame_route=True,
    inverse_recurrences_and_cutoffs_are_coefficientwise_monotone_in_reciprocal=True,
    full_moving_normal_unprojected_pressure_source_and_global_memory_retained=True,
    original_primitive_validity_physical_domain_and_current_source_obligations_retained=True,
    new_analytical_assumptions_or_selector_thresholds_added=False,
    actual_primitive_enclosures_unknown_historical_band_or_generated_source_bounds_certified=False,
    new_Mz_Db_or_CE_substituted_into_old_length=False,
    roots_fields_carriers_tilt_support_width_or_historical_evidence_changed=False,
    unconditional_pulse_complete_residual_infinite_field_or_blowup_proved=False)
REMAINING=[
    'retain validity of the ORIGINAL full background/frame primitive bounds and same physical accepted-tail domain; native leading first jets alone do not prove Db or every phase jet',
    'establish current generated complete harmonic full-path source/support and weighted target/freezing/root/cone inequalities using the already derived normal and tangent implications',
    'close full covariance and fresh axial/full-memory pressure means plus every viscosity/cutoff/interaction residual',
    'verify infinite field localization forcing and complete velocity growth before unconditional pulse or blow-up claims']
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
    paths += [HERE.parent/p for p in ('core_normal_length_closure/budgets.py',
        'core_sourced_normal_transfer/budgets.py','core_tangent_propagator_closure/budgets.py',
        'frame_column_foundation/frame_budget.py','core_first_sourced_inverse_budget/budgets.py',
        'core_leading_direction_enclosure/bounds.py','actual_band_pulse_bounds/README.md')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {str(k):canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_reference_floor_and_inverse_controls=bool(r['checks']) and all(r['checks'].values()),
            eight_historical_records_bound=len(bound)==8,
            all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Reference-scaled original normal floor, inverse transfer, frontier or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
