#!/usr/bin/env python3
"""Freeze same-gap reweighting and fresh post-sourced mean SIZE closure."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_sourced_covariance_moment_join','core_signed_prefix_source_budget','core_initial_haar_row_budget','core_dynamic_mean_recheck')
STATUS='same-gap-higher-pole-operators-and-fresh-post-sourced-full-pressure-mean-size-closure'
FLAGS={
    'higher_pole_operator_multipliers_derived_from_retained_positive_profile_recipes':True,
    'actual_original_gap_width_and_closed_reserved_hulls_retained':True,
    'finite_shift_is_estimate_argument_only_without_new_physical_band_gate':True,
    'old_final_means_complete_changed_wave_and_fresh_fast_slow_maps_composed':True,
    'complete_pressure_and_three_scalars_recomputed_prefast_postfast_and_postslow':True,
    'next_source_interface_rejects_inherited_single_frequency_grammar':True,
    'new_analytical_assumptions_or_acceptance_gates_added':False,
    'arbitrary_old_operator_norm_implies_new_domain_coverage':False,
    'actual_original_input_enclosures_or_profile_provenance_certified':False,
    'homogeneous_frame_or_original_root_identification_closed':False,
    'complete_residual_improvement_unconditional_pulses_or_blowup_proved':False,
    'historical_evidence_profiles_supports_widths_or_bands_changed':False,
}
REMAINING=[
    'validate original primitive old-wave source common-chart coefficient/support and recipe-defined profile bounds across every actual physical scale; no manufactured fixture certifies these inputs',
    'derive exact next nonzero regular source for the changed +/-1,+/-2 wave and fresh means including all output harmonics through +/-4 and full pressure angular factors',
    'establish the complete improved residual inequality including viscosity actual cutoff commutators tails all interactions and fresh zero-mode identities on the original band',
    'retain and where needed identify original homogeneous primary frames roots actual A0/A1 velocities and moving-pressure compensator products',
]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks


def sources():
    paths=[HERE/n for n in ('README.md','budgets.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/n/'budgets.py' for n in ('core_sourced_covariance_moment_join','core_signed_prefix_source_budget','core_initial_haar_row_budget','core_signed_target_edge_budget','core_tangent_propagator_closure')]
    paths += [HERE.parent/'core_primary_aggregate_atlas/proof_checks.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_kernel_shift_and_fresh_mean_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   four_historical_records_bound=len(bound)==4,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Reweighted mean size derivation, retained recipe provenance or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} inherited source checks")


if __name__=='__main__':main()
