#!/usr/bin/env python3
"""Freeze canonical principal Haar-column jets, geometry and original-root norm join."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_canonical_entrance_jet_closure','gaussian_mass_foundation',
    'auxiliary_geometry_audit','core_phase_coefficient_budget','core_primary_root_budget',
    'canonical_primary_foundation','frame_column_foundation','core_primary_aggregate_atlas')
STATUS='conditional-canonical-Haar-column-jet-closure'
FLAGS=dict(original_full_principal_column_integrals_and_six_tensor_entries_retained=True,
    original_angular_half_transverse_mass_Jacobian_and_physical_Haar_scale_retained=True,
    moving_translation_derivatives_cancel_only_in_complete_normalized_Haar_integral=True,
    finite_common_cover_preserves_normalized_Haar_without_deck_count=True,
    raw_column_jets_derived_from_same_canonical_pulse_and_Gaussian_square_bound=True,
    actual_fixed_positive_tilt_not_inferred_from_an_upper_bound=True,
    root_H_and_W_norm_inputs_jointly_derived_without_rebuilding_physical_roots=True,
    full_forced_pressure_memory_axial_mean_and_reserved_support_obligations_retained=True,
    canonical_and_inherited_families_or_principal_and_full_curl_covariances_identified=False,
    actual_full_background_suprema_weighted_target_or_full_curl_operator_input_certified=False,
    original_band_length_width_carriers_profiles_fields_or_historical_evidence_changed=False,
    unconditional_pulse_complete_residual_infinite_field_or_blowup_proved=False)
REMAINING=[
    'derive full curl-operator D_j from actual normal reciprocal, cylindrical/graph/potential coefficients and all coefficient derivatives, or state the exact blocking inequality',
    'retain full actual background/phase jet validity, weighted target/freezing/cone/value hypotheses and original accepted physical domains',
    'join all complete-wave covariance means with axial velocity, full pressure primitives/moments and reserved supports',
    'complete every viscosity/cutoff/interaction residual and infinite realization before unconditional pulse, complete residual improvement or blow-up claims']
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
    paths += [HERE.parent/p for p in ('core_canonical_entrance_jet_closure/bounds.py',
        'gaussian_mass_foundation/gaussian_bounds.py','auxiliary_geometry_audit/geometry.py',
        'core_phase_coefficient_budget/budgets.py','core_primary_root_budget/budgets.py',
        'actual_band_pulse_bounds/README.md','canonical_initial_mean_join/README.md')]
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {str(k):canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_column_geometry_and_joint_root_controls=bool(r['checks']) and all(r['checks'].values()),
            eight_historical_records_bound=len(bound)==8,
            all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Canonical column integral, Haar geometry, joint root norm recipe, scope or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
