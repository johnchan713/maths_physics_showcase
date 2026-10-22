#!/usr/bin/env python3
"""Freeze generated complete covariance, pressure moments and fresh pole domains."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_tangent_propagator_closure','core_signed_prefix_source_budget','core_dynamic_mean_recheck','first_sourced_foundation')
STATUS='generated-sourced-complete-covariance-full-pressure-moments-and-next-weighted-operator-domain-correction'
FLAGS={
    'new_covariance_bound_generated_from_before_source_derived_inverse_and_same_complete_cut_curl':True,
    'all_signed_modes_quadratures_axial_entries_and_second_harmonic_square_retained':True,
    'three_changed_defects_and_whole_moving_axial_pressure_moment_bounded':True,
    'full_covariance_support_hull_not_reserved_mean_hull_used':True,
    'next_weighted_operator_domains_6p_and_12p_kept_distinct_from_old_2p_and_4p':True,
    'new_analytical_assumptions_or_acceptance_gates_added':False,
    'actual_primitive_old_wave_or_before_source_enclosures_certified':False,
    'actual_new_weighted_operator_constants_and_post_wave_fresh_mean_fields_certified':False,
    'homogeneous_primary_frame_or_original_root_identification_closed':False,
    'complete_residual_improvement_unconditional_pulses_or_blowup_proved':False,
    'historical_profiles_supports_widths_bands_or_evidence_changed':False,
}
REMAINING=[
    'retain validity of original primitive old-wave source coefficient and full-path support bounds on the same common chart; generated covariance is no longer an independent norm input on this route',
    'reevaluate SAME weighted operators at 6p and 12p (or prove sharper actual poles), then apply fresh fast and slow maps and recompute complete pressure and all three scalars after each update',
    'retain old means actual A0/A1 velocities full Pi1/Pi2 memory moving compensator products and homogeneous primary/root identification where required',
    'join full target and pointwise residual estimates across all physical scales including viscosity cutoff and interaction terms; scalar moment bounds alone do not prove residual improvement',
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
    paths += [HERE.parent/'core_tangent_propagator_closure/budgets.py',HERE.parent/'core_signed_prefix_source_budget/budgets.py',
              HERE.parent/'core_signed_target_edge_budget/budgets.py',HERE.parent/'core_primary_aggregate_atlas/proof_checks.py',
              HERE.parent/'complete_covariance_initial_means/means.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_tensor_curl_pressure_pole_and_generated_chain_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   four_historical_records_bound=len(bound)==4,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Generated covariance pressure moments weighted proof domains or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} inherited source checks")


if __name__=='__main__':main()
