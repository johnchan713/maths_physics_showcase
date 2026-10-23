#!/usr/bin/env python3
"""Freeze regional leading enclosures and the unchanged-selector obstruction."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_native_axial_transport_floor','higher_jet_budget',
    'global_axial_source_envelope','global_swirl_envelope','collar_operator_envelope',
    'global_axial_angular_bound','leading_frame_foundation','background_jet_foundation',
    'core_primary_tail_recheck','core_physical_representative_closure')
STATUS='regional-leading-first-jets-joined-with-same-selector-margin-open'
FLAGS=dict(actual_defining_map_leading_direction_bounds_analytically_assembled=True,
    first_compact_logarithmic_jets_have_no_N_loss=True,
    second_radial_jets_still_have_N_loss=True,
    whole_physical_hull_covered_by_original_regional_patches=True,
    regular_swirl_half_power_and_one_sided_angular_jets_retained=True,
    same_complete_first_Euler_diagonal_and_primary_tail_slow_gate_identified=True,
    full_pressure_memory_radial_cutoff_curl_and_same_roots_retained=True,
    actual_giant_parameters_or_same_Lreq_margin_numerically_verified=False,
    unknown_historical_first_band_identified_with_primary_tail=False,
    new_selector_Lreq_tilt_width_frequency_or_support_selected=False,
    inherited_continuum_premises_independently_verified_unconditional=False,
    unconditional_pulse_complete_residual_infinite_field_or_blowup_proved=False,
    historical_evidence_changed=False)
REMAINING=[
    'evaluate the explicit same-profile leading/reference/geometry expressions in 4[A/Sfloor+B/Sfloor^2+C/Sfloor^3]^2<=beta for the EXISTING Lreq',
    'verify same-band target freezing, original root/cone and generated complete harmonic path/source norms',
    'join complete-wave covariance axial/pressure mean updates and every viscosity/cutoff/interaction residual',
    'verify infinite field localization forcing and full velocity growth before any blow-up claim']
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
    paths += [HERE.parent/'core_native_axial_transport_floor/bounds.py',
        HERE.parent/'leading_frame_foundation/frame_floors.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {str(k):canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_leading_enclosure_controls=bool(r['checks']) and all(r['checks'].values()),
            ten_historical_records_bound=len(bound)==10,
            all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Regional leading enclosure, unchanged-selector scope or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
