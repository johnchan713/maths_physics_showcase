#!/usr/bin/env python3
"""Freeze localized canonical stress assembly and its stage-specific row join."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_haar_row_stress_join','core_fast_row_gain_join','canonical_primary_foundation',
              'auxiliary_geometry_audit','actual_band_pulse_bounds','current_residual_ledger_audit')
STATUS='conditional-localized-canonical-principal-stress-and-complete-regular-postwave-tangential-D'
FLAGS={
    'both_pulse_cutoffs_angular_half_and_Haar_measure_inside_same_columns':True,
    'canonical_positive_roots_and_physical_squared_partition_match_same_stress':True,
    'three_term_defect_identity_retains_root_scale_and_cover_errors':True,
    'exact_partition_identity_differentiated_after_entire_assembly':True,
    'both_signs_counted_inside_each_box_matrix_not_twice_in_cover':True,
    'full_regular_postwave_tangential_mean_rows_D_derived_conditionally':True,
    'precise_postwave_stage_full_pressure_and_high_jet_A1_recipe_retained':True,
    'spatial_cutoffs_not_reclassified_as_Gaussian_flat_tails':True,
    'final_Haar_C_and_radial_Haar_B_plus_C_preserved':True,
    'new_analytical_premise_or_band_gate_added':False,
    'independent_actual_column_root_profile_or_support_inputs_certified':False,
    'two_different_homogeneous_families_identified_by_closeness':False,
    'postwave_rows_identified_with_final_fresh_rows':False,
    'full_pointwise_radial_A0_or_complete_physical_tail_bound_proved':False,
    'complete_residual_unconditional_pulses_infinite_field_or_blowup_proved':False,
    'historical_frozen_evidence_fields_or_width_changed':False,
}
REMAINING=[
    'trace final fresh feedback through full tangential rows and bound the pointwise radial A0 remainder without promoting its Haar estimate',
    'join every exact Gaussian/spatial/scale cutoff and interaction tail with complete wave pressure and viscosity to the physical residual',
    'independently validate original actual background/cone/column/root/profile/support/endpoint contracts including the very high finite A1 jet inputs',
    'prove uniform actual pulses, forcing regularity, infinite realization and complete velocity growth before any blow-up claim',
]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    paths=[HERE/n for n in ('README.md','assembly.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_haar_row_stress_join/budgets.py',HERE.parent/'core_fast_row_gain_join/budgets.py',
              HERE.parent/'core_reweighted_fresh_mean_closure/budgets.py',HERE.parent/'core_primary_root_budget/README.md',
              HERE.parent/'canonical_initial_mean_join/README.md',HERE.parent/'core_initial_haar_row_budget/README.md']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_localization_and_stage_composition_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   six_historical_records_bound=len(bound)==6,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Localized stress, stage-specific row join, historical provenance or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")


if __name__=='__main__':main()
