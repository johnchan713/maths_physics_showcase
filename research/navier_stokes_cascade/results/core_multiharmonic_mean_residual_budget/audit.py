#!/usr/bin/env python3
"""Freeze multiharmonic regular source and the actual mean-feedback gain gap."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_reweighted_fresh_mean_closure','core_sourced_covariance_moment_join','core_signed_prefix_source_budget',
              'regular_source_foundation','current_residual_ledger_audit','finite_induction_foundation')
STATUS='quantitative-multiharmonic-current-regular-source-and-sharp-coarse-mean-feedback-gain-gap'
FLAGS={
    'quantitative_next_source_uses_all_signed_first_and_second_input_harmonics':True,
    'input_ratio_two_and_principal_damping_m_squared_four_retained':True,
    'full_uncut_source_normal_pressure_and_actual_complete_cut_velocities_retained':True,
    'exact_wave_cancellation_and_both_fresh_mean_wave_feedback_orientations_retained':True,
    'zero_mode_full_pressure_and_all_three_scalars_reconstructed':True,
    'explicit_coarse_order_B_source_and_nine_pole_inverse_interface_derived':True,
    'exact_additive_tail_ledger_kept_outside_every_inverse_input':True,
    'actual_mean_feedback_order_D_is_precise_remaining_inequality':True,
    'class_level_phase_advection_witness_shows_coarse_size_alone_insufficient':True,
    'new_analytical_assumptions_or_accepted_band_gates_added':False,
    'fresh_actual_mean_feedback_improved_order_certified':False,
    'witness_is_the_prescribed_actual_fast_or_slow_map_output':False,
    'historical_conditional_cycle_or_actual_maps_refuted':False,
    'all_spatial_localization_errors_are_Gaussian_tails':False,
    'actual_original_input_enclosures_or_full_support_contract_certified':False,
    'complete_residual_improvement_unconditional_pulses_infinite_field_or_blowup_proved':False,
    'historical_evidence_width_profiles_carriers_supports_or_band_changed':False,
}
REMAINING=[
    'prove actual fresh mean/current-wave feedback inequality (8) at D, or establish sufficient actual post-wave zero-Haar fast row D and fresh postfast scalar C bounds using exact current targets; coarse mean size alone does not imply them',
    'retain complete zero-mode pressure and moment identities with actual A0/A1 fields and whole moving-compensator products when proving any complete residual gain',
    'validate original actual primitive old-wave source coefficient/profile recipe common-chart and closed support inputs uniformly across every physical scale and identify original primary/root frames where needed',
    'join any eventual regular gain to unchanged exact additive tails, viscosity and spatial/diagonal cutoff estimates and the infinite-field/forcing obligations without forwarding flat tails to an inverse',
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
    paths += [HERE.parent/n/'budgets.py' for n in ('core_reweighted_fresh_mean_closure','core_sourced_covariance_moment_join','core_signed_prefix_source_budget')]
    paths += [HERE.parent/'core_reweighted_fresh_mean_closure/controls.py',HERE.parent/'regular_source_foundation/source_checks.py',
              HERE.parent/'current_residual_ledger_audit/residual_ledger.py',HERE.parent/'finite_induction_foundation/finite_checks.py',
              HERE.parent/'first_sourced_foundation/sourced_checks.py',HERE.parent/'complete_covariance_initial_means/means.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}


def generate():
    r=review();bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_harmonic_witness_polynomial_and_composed_controls_pass=bool(r['checks']) and all(r['checks'].values()),
                   six_historical_records_bound=len(bound)==6,all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Multiharmonic source, mean-feedback gap, provenance or scientific scope changed')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} inherited source checks")


if __name__=='__main__':main()
