#!/usr/bin/env python3
"""Freeze final regular angular-zero closure without promoting actual inputs."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_localized_stress_assembly','core_haar_row_stress_join','core_fast_row_gain_join',
              'core_fast_inverse_budget','core_actual_scalar_gain_join','core_reweighted_fresh_mean_closure','compact_mean_foundation')
STATUS='conditional-final-regular-angular-zero-rows-C-and-compensated-pointwise-A0-bound'
FLAGS={
    'fresh_fast_cancels_old_Q1_rows_and_retains_exact_A1_error':True,
    'complete_pressure_viscosity_and_actual_mean_differences_stronger_than_C':True,
    'fresh_A1_uses_two_integrations_and_retains_nested_old_high_jets':True,
    'final_tangential_Q1_three_halves_and_Haar_C_joined':True,
    'unit_mass_compensation_kills_integrated_A0_zero_mode_exactly':True,
    'pointwise_A0_is_bounded_not_deleted_by_Haar_cancellation':True,
    'current_full_graph_radial_source_and_original_e0_kernel_profiles_retained':True,
    'final_all_three_regular_mean_rows_C_and_combined_regular_D_derived':True,
    'proper_poles_width_gap_and_exact_additive_tail_ledger_unchanged':True,
    'new_analytical_assumption_or_starting_band_gate_added':False,
    'pointwise_pressure_field_promoted_to_C':False,
    'actual_high_jet_polynomials_or_full_profile_suprema_evaluated_and_verified':False,
    'complete_physical_tail_or_full_residual_bound_proved':False,
    'unconditional_actual_pulses_infinite_realization_or_blowup_proved':False,
    'historical_frozen_evidence_fields_or_width_changed':False,
}
REMAINING=[
    'join every exact additive base, Gaussian, spatial and scale cutoff tail with complete wave pressure, viscosity and interactions after physical losses',
    'independently validate actual background/cone/root/profile and support/endpoint contracts including current A0 source jets through 100005 and nested A1 original mean jets through 200012 at output zero',
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
    paths=[HERE/n for n in ('README.md','budgets.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/n/'budgets.py' for n in ('core_haar_row_stress_join','core_fast_row_gain_join','core_actual_scalar_gain_join','core_reweighted_fresh_mean_closure','core_signed_prefix_source_budget')]
    paths += [HERE.parent/'core_localized_stress_assembly/assembly.py',HERE.parent/'core_fast_inverse_budget/proof_checks.py',
              HERE.parent/'core_reweighted_fresh_mean_closure/controls.py',HERE.parent/'core_multiharmonic_mean_residual_budget/README.md']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_compact_pressure_fresh_difference_and_composition_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   seven_historical_records_bound=len(bound)==7,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Final regular mean closure, pointwise A0, high-jet limitations, provenance or scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")


if __name__=='__main__':main()
