#!/usr/bin/env python3
"""Freeze fixed-profile target, flat-log, cutoff composition and same-root joins."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_complete_curl_operator_budget','weighted_background_edges','weighted_stress_foundation',
    'core_fixed_profile_tame','auxiliary_geometry_audit','core_localized_stress_assembly',
    'core_primary_root_budget','core_initial_haar_row_budget')
STATUS='conditional-fixed-profile-primary-target-and-cutoff-budget'
FLAGS=dict(original_inner_outer_factorizations_give_smooth_cubic_pole_numerator=True,
    actual_leading_target_and_scale_prefactor_retained_with_sharper_n_plus_three_pole=True,
    nonlinear_profile_chart_and_all_positive_flat_log_jets_retained=True,
    original_flat_width_square_and_outer_constant_four_retained=True,
    fixed_cutoff_profile_seminorms_and_original_S_cubed_mesh_supplied=True,
    relative_temporal_cutoff_uses_exact_ci_Ls_equals_two_r0=True,
    slow_partition_included_once_in_wave_norm_and_kept_outside_H=True,
    every_primary_norm_slot_joined_on_same_original_value_and_cone_domain=True,
    all_pressure_memory_axial_means_reserved_supports_and_residual_terms_retained=True,
    actual_giant_profile_seminorms_full_background_or_value_gate_numerically_certified=False,
    original_width_radius_carriers_profiles_roots_or_historical_evidence_changed=False,
    unconditional_pulse_full_physical_residual_infinite_field_or_blowup_proved=False)
REMAINING=[
    'substitute explicit primary inputs into existing scalar/source and complete covariance/mean joins on the SAME weighted physical hull, retaining their linked prerequisites',
    'retain validity of original full background/phase jet polynomials, fixed defining-map/profile seminorms and original physical value/cone domain; no numerical giant parameter certification follows',
    'retain all six covariance entries, fresh axial velocity, Pi1/Pi2 pressure memory, moving moments and reserved supports in every mean update',
    'identify and close the exact remaining viscosity/cutoff/interaction full residual inequality and infinite realization before unconditional pulse or blow-up claims']
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
    paths += [HERE.parent/p for p in ('core_complete_curl_operator_budget/bounds.py',
        'core_canonical_haar_column_jet_closure/bounds.py','core_canonical_entrance_jet_closure/bounds.py',
        'core_fixed_profile_tame/proof_checks.py','core_reserved_geometry_recheck/proof_checks.py',
        'auxiliary_geometry_audit/geometry.py','core_primary_root_budget/budgets.py',
        'primary_covariance_foundation/README.md','core_native_chart_cover/README.md')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {str(k):canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_target_chart_cutoff_and_join_controls=bool(r['checks']) and all(r['checks'].values()),
            eight_historical_records_bound=len(bound)==8,
            all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Fixed-profile target/loghalf/cutoff, chart, original-root join or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
