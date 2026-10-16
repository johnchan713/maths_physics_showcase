#!/usr/bin/env python3
"""First sourced core correction with exact excluded Gaussian cutoff tails."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from sources import load,cutoff_review,forced_normal_review,tail_review,covariance_pressure_review,exponent_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='first-sourced-wave-core-gain-with-exact-retained-cutoff-tails'
FLAGS=dict(current_regular_source_class_support_and_endpoint_hypotheses_joined=True,
           first_sourced_inverse_inputs_and_all_fixed_physical_jets_joined=True,
           exact_full_potential_pressure_cutoff_commutator_retained=True,
           homogeneous_flat_source_nonflat_response_obstruction_derived=True,
           all_six_complex_covariances_and_full_axial_pressure_recomputed=True,
           mean_residual_and_three_defects_regenerated_not_preserved=True,
           one_finite_regular_normalized_residual_order_gain_joined=True,
           both_base_pressure_tails_old_Haar_means_and_actual_flat_velocities_retained=True,
           relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas=True,
           manufactured_diagnostics_not_actual_giant_band_evaluation=True,
           forward_inverse_asserted_to_preserve_physical_flatness=False,
           zero_entry_response_asserted_zero_at_exit=False,
           designated_flat_tails_forward_inverted_or_erased=False,
           initial_real_quadrature_shortcut_used_for_complex_update=False,
           actual_giant_profile_or_first_band_numerically_evaluated=False,
           every_later_complex_wave_source_and_full_cycle_verified=False,
           unconditional_complete_residual_improvement_verified=False,
           complete_infinite_potential_realization_verified=False,
           blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['apply_generic_complete_wave_signed_inverse_with_current_remainder_rho_1_minus_3kappa',
           'relax_regenerated_full_means_and_defects_recomputing_pressure_after_each_map',
           'verify_all_later_sources_supports_and_complete_finite_cycles_with_exact_tail_ledger',
           'establish_infinite_realization_forcing_and_energy_before_any_blowup_claim']
DEPENDENCIES=('post_signed_mean_relaxation','real_signed_covariance_update','physical_pulse_transfer',
              'actual_band_pulse_bounds','geometric_pulse_replay','residual_cycle_closure',
              'compact_mean_corrections','physical_wave_interface')
REVIEWERS=dict(full_viscous_potential_pressure_cutoff=cutoff_review,
               constrained_sourced_pressure_and_endpoint=forced_normal_review,
               exact_forward_inverse_flat_tail_obstruction=tail_review,
               general_complex_covariance_and_axial_pressure=covariance_pressure_review,
               regular_source_output_and_physical_ledger=exponent_review)
EXTRA_GATES={'current_post_signed_parent_replayed','frozen_provenance_preserved'}


def protocol():
    return dict(parent_commit='6eeb414542ecf473826a6bc9174fc050fd73160b',
        paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_equations=['6.19-6.32','7.39-7.42','8.12-8.27','9.1-9.10'],
        same_Md64_h_N_positive_lambda_width_and_global_geometry=True,
        kappa='1/100000',B='1-3*kappa',C='3/2-3*kappa',D='3/2-6*kappa',mu='1-kappa',
        current_harmonics=[-2,-1,1,2],updated_residual_maximum_multiplier=4,
        exact_source_split='G_total=G_core+T_old; invert only G_core',
        operation_order='split exact tails; sourced inverse; cut potential AND pressure; recompute full covariance and mean pressure',
        source_smooth_zero_extension_not_assumed_for_uncut_output_exit=True,
        designated_Gaussian_flat_tails_never_forward_inverted=True,
        whole_moving_pressure_bump_product_differentiated=True,
        all_six_covariances_both_base_tails_and_actual_A0_A1_retained=True,
        updated_fixed_primary_remainder_order='1-3*kappa',
        next_signed_inverse_requires_generic_complex_covariance=True,
        fixture='manufactured compact physical jets; exact rational moving normal and complex quadratures',
        diagnostic_step_extension='C2; analytical obstruction uses genuine C-infinity step',
        physical_sample_h='1/200',enlarged_chart_domain='4Q<qstar',
        finite_stage_order_gain_not_uniform_iteration_or_fixed_scale_contraction=True,
        diagnostics_not_continuum_proofs=True)


@lru_cache(None)
def parent_module():
    p=HERE.parent/'post_signed_mean_relaxation';sys.path.insert(0,str(p))
    try: return load('first_source_parent_post_signed',p/'audit.py')
    finally: sys.path.pop(0)


@lru_cache(None)
def parent_review():
    parent=parent_module();record=parent.generate();parent.validate(record)
    frozen=json.loads((HERE.parent/'post_signed_mean_relaxation/evidence.json').read_text())
    parent.compare_records(record,frozen)
    inherited=dict(record['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['post_signed_mean_relaxation:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return record,inherited


def source_hashes():
    paths=[HERE/n for n in ('README.md','sources.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol(): raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review();reviews={n:f() for n,f in REVIEWERS.items()}
    gates={n:bool(r['checks']) and all(r['checks'].values()) for n,r in reviews.items()}
    gates['current_post_signed_parent_replayed']=all(parent['gates'].values()) and all(parent['flags'][n] for n in
        ('current_post_signed_fast_and_slow_mean_inputs_joined',
         'full_angular_mean_tangential_residual_order_3_minus_6kappa',
         'complete_axial_pressure_moment_and_both_base_tails_retained'))
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'): return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:
        raise ValueError('Changed construction, exact tail split or restricted scope')
    if r.get('remaining_obligations')!=REMAINING: raise ValueError('Outstanding obligation changed')
    if set(r.get('reviews',{}))!=set(REVIEWERS) or set(r.get('gates',{}))!=set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(r['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in r['reviews'].values()):
        raise ValueError('First sourced wave check failed')
    if len(r.get('inherited_provenance',{}))!=640 or not all(r['inherited_provenance'].values()):
        raise ValueError('Inherited provenance changed')
    if r.get('source_hashes')!=source_hashes(): raise ValueError('Current sources changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')
    for n,f in REVIEWERS.items(): compare_records(f(),r['reviews'][n],n)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args(argv)
    if args.output.exists(): raise SystemExit('Output exists; choose a fresh path')
    r=generate();validate(r)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text());validate(frozen);compare_records(r,frozen)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(STATUS)
    print(f"{len(r['gates'])} gates pass; {len(r['inherited_provenance'])} inherited hashes agree.")
    print('One finite regular-source step with exact flat tails; later cycles/full construction remain open.')


if __name__=='__main__': main()
