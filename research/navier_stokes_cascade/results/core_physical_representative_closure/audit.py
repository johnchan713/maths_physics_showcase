#!/usr/bin/env python3
"""Freeze original representative membership and the physical norm-only floor."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_native_chart_cover','core_tangential_normal_closure','core_borel_prescription_closure','auxiliary_geometry_audit','leading_frame_foundation')
STATUS='original-physical-representative-membership-and-norm-only-tangential-floor-closure'
FLAGS={
    'original_representatives_are_physical_active_points_by_defining_rule':True,
    'existing_physical_doubled_box_norm_hull_applies_from_original_representative':True,
    'reference_growth_shear_gaps_needed_at_active_representative_only_for_this_floor':True,
    'independent_native_value_neighborhood_and_negative_time_rep_requirements_removed_for_floor':True,
    'complete_axial_transport_pressure_memory_and_radial_cutoff_premises_retained':True,
    'new_analytical_assumptions_added':False,
    'global_negative_time_extension_or_full_cone_neighborhood_closed':False,
    'pdf_binary_hash_independently_recomputed':False,
    'actual_profile_input_constants_or_unchanged_band_floor_numerically_verified':False,
    'historical_representatives_scales_supports_or_evidence_changed':False,
    'complete_residual_unconditional_pulses_or_blowup_proved':False,
}
REMAINING=[
    'substitute actual leading first tuple norms on original physical norm hull and closed-active-annulus reference gap bounds into unchanged-band normal and slow inequalities',
    'retain original radial clearance and q-domain gates; stricter supplied rational root enclosures need their own clearance check',
    'keep native value neighborhood growth/cone estimates and global smooth extension obligations where OTHER full-wave estimates require them',
    'discharge remaining sourced propagator supports full covariance fresh pressure/axial means and complete viscous cutoff interaction residual estimates',
]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    paths=[HERE/n for n in ('README.md','budgets.py','controls.py','audit.py','test_audit.py','protocol.json','source_contract.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_native_chart_cover/budgets.py',HERE.parent/'core_tangential_normal_closure/budgets.py',HERE.parent/'background_jet_foundation/protocol.json']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),source_contract=json.loads((HERE/'source_contract.json').read_text()),flags=deepcopy(FLAGS),review=r,
                gates=dict(all_exact_representative_norm_and_floor_controls_pass=bool(r['checks']) and all(r['checks'].values()),five_historical_records_bound=len(bound)==5,
                           all_declared_historical_sources_match=bool(checks) and all(checks.values())),dependency_evidence_sha256=bound,
                inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Original physical representative floor-only implication or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} inherited source checks")


if __name__=='__main__':main()
