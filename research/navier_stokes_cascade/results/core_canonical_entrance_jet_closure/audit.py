#!/usr/bin/env python3
"""Freeze the canonical entrance, homogeneous jet closure and unchanged root join."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_reference_scaled_normal_floor','core_tangent_propagator_closure',
    'core_first_sourced_inverse_budget','core_phase_coefficient_budget',
    'core_primary_root_budget','canonical_primary_foundation','frame_column_foundation',
    'core_initial_haar_row_budget')
STATUS='conditional-canonical-entrance-and-homogeneous-jet-closure'
FLAGS=dict(original_nonzero_canonical_entrance_and_reference_mixing_retained=True,
    tangential_floor_uses_original_full_reference_error_proof=True,
    Gaussian_initial_derivatives_and_moving_normal_pressure_retained=True,
    original_accepted_tail_tangent_bound_supplies_all_mixed_homogeneous_jets=True,
    independent_homogeneous_wave_jet_callback_replaced_by_lower_input_recipe=True,
    root_join_changes_norm_slot_only_and_keeps_cutoff_derivatives=True,
    fixed_entrance_label_and_full_actual_phase_jet_hypotheses_retained=True,
    complete_forced_pressure_source_and_global_memory_obligations_retained=True,
    initial_weighted_row_and_generated_source_conditional_joins_not_reopened=True,
    canonical_and_inherited_physical_primaries_identified=False,
    new_selector_roots_fields_carriers_tilt_or_support_width_changed=False,
    actual_background_suprema_Haar_column_jets_or_all_generated_sources_certified=False,
    unconditional_pulse_complete_residual_infinite_field_or_blowup_proved=False)
REMAINING=[
    'derive canonical Haar-column H_j jets from the homogeneous recipe and full covariance integrals with moving maps, cutoffs, Gaussian and edge normalization, or state the exact blocking inequality',
    'retain validity of full actual phase/background jets and original common physical accepted domains; leading first jets alone are insufficient',
    'supply every lower input in the existing generated-source/weighted-target and finite support joins, retaining global pressure/moment and reserved-mean supports',
    'close fresh complete-wave covariance means including axial and full-memory pressure, then every viscosity/cutoff/interaction residual before any unconditional pulse or blow-up conclusion']
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
    paths += [HERE.parent/p for p in ('core_reference_scaled_normal_floor/bounds.py',
        'core_sourced_normal_transfer/budgets.py','core_tangent_propagator_closure/budgets.py',
        'core_first_sourced_inverse_budget/budgets.py','core_phase_coefficient_budget/budgets.py',
        'core_primary_root_budget/budgets.py','core_signed_prefix_source_budget/README.md')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {str(k):canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_canonical_entrance_and_homogeneous_jet_controls=bool(r['checks']) and all(r['checks'].values()),
            eight_historical_records_bound=len(bound)==8,
            all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Canonical entrance, homogeneous jet recurrence, root join, scope or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
