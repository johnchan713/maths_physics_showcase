#!/usr/bin/env python3
"""Replay physical primary pulse transfer; nonlinear source estimates stay open."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from transfer import load,exact_review,threshold_review,graph_review,tail_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='actual-background-to-physical-primary-pulse-transfer-audited'
FLAGS=dict(actual_background_inputs_to_primary_pulse_theorem_discharged=True,
           common_all_later_band_primary_pulse_values_and_fixed_jets_derived=True,
           complete_primary_curl_velocity_and_wave_pressure_physically_transferred=True,
           exact_closed_support_and_one_sided_boundary_conditions_retained=True,
           sourced_inverse_requires_input_class_and_smooth_endpoint_extension=True,
           first_two_mean_pressure_tails_and_complete_covariance_obligation_retained=True,
           relies_on_replayed_background_geometry_and_pulse_lemmas=True,
           actual_giant_profile_cutoff_scales_or_first_band_numerically_evaluated=False,
           all_nonlinear_generated_source_class_estimates_verified=False,
           unconditional_complete_nonlinear_pulse_hierarchy_verified=False,
           complete_covariance_mean_update_verified=False,
           unconditional_complete_residual_improvement_verified=False,
           blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['audit_complete_wave_covariance_mean_updates_with_axial_velocity_and_pressure_tails',
           'verify_class_and_support_estimates_for_every_nonlinear_generated_source',
           'prove_complete_residual_gain_with_viscosity_cutoffs_and_all_interactions']
DEPENDENCIES=('weighted_background_edges','geometric_pulse_replay','actual_band_pulse_bounds',
              'auxiliary_geometry_audit','physical_wave_interface','background_jet_foundation',
              'coefficient_support_audit','compact_mean_corrections')
REVIEWERS=dict(exact_graph_and_physical_power_ledger=exact_review,
               common_band_domain_and_enlargement=threshold_review,
               independent_physical_graph_and_phase_derivatives=graph_review,
               physical_flat_tails_and_closed_supports=tail_review)


def protocol():
    return dict(parent_commit='7c325b42d7e517d60a9208bd562012a4fb357f1b',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                same_Md64_h_N_width_and_global_geometry=True,radius='1/6057072',
                spatial_harmonic_loss='1/2+h/2',temporal_harmonic_loss='1+3h/2',
                common_potential_loss='2A+(m+1)*(1+3h/2)',
                enlarged_background_margin='4Q<qstar',working_domain='q<=2^(-ell0-2)',
                all_later_bands_not_a_finite_band_scan=True,
                actual_first_band_and_coefficient_constants_numerically_evaluated=False,
                fixed_labels_carriers_covers_Q_epsilon_S_during_differentiation=True,
                derivative_graph='radial graph cancels in phase but not amplitudes; fast auxiliary time retained',
                physical_rescalings=['Q^-A velocity','Q^-2A pressure','Q^(1/2-A) fixed vector potential'],
                primary_normalized_exponents=['1/2 velocity','1 pressure','1 potential'],
                source_inverse_has_explicit_unverified_input_class_obligation=True,
                complete_covariance_means_and_full_residual_not_inferred=True,
                diagnostic_graph_bands=[160,1024],diagnostic_h='1/200',
                manufactured_derivative_lift_translation_not_actual_map_change=True,
                diagnostics_are_not_continuum_proofs=True)


@lru_cache(None)
def parent_module():
    p=HERE.parent/'weighted_background_edges'
    sys.path.insert(0,str(p))
    try:
        return load('physical_transfer_parent_weighted_audit',p/'audit.py')
    finally:
        sys.path.pop(0)


def parent_review():
    module=parent_module();record=module.generate();module.validate(record)
    frozen=json.loads((HERE.parent/'weighted_background_edges/evidence.json').read_text())
    module.compare_records(record,frozen)
    inherited=dict(record['inherited_provenance'])
    for path,expected in frozen['source_hashes'].items():
        inherited['weighted_background_edges:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==expected
    return record,inherited


def source_hashes():
    paths=[HERE/n for n in ('README.md','transfer.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():
        raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review()
    reviews={k:f() for k,f in REVIEWERS.items()}
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews.items()}
    gates['audited_background_inputs_replayed']=all(parent['gates'].values()) and all(parent['flags'][k] for k in
        ('leading_all_fixed_order_weighted_edge_certificate_joined','actual_background_relative_stress_comparison_derived',
         'actual_background_growth_and_oriented_cone_on_common_small_q_domain'))
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
                inherited_provenance=inherited,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):
    return parent_module().compare_records(a,b,path)


def validate(record):
    if record.get('status')!=STATUS or record.get('protocol')!=protocol() or record.get('flags')!=FLAGS:
        raise ValueError('Changed construction or scope')
    if record.get('remaining_obligations')!=REMAINING:
        raise ValueError('Outstanding obligation removed')
    if set(record.get('reviews',{}))!=set(REVIEWERS) or set(record.get('gates',{}))!=set(REVIEWERS)|{'audited_background_inputs_replayed','frozen_provenance_preserved'}:
        raise ValueError('Review or gate removed')
    if not all(record['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in record['reviews'].values()):
        raise ValueError('Physical transfer check failed')
    if len(record.get('inherited_provenance',{}))!=580 or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited provenance changed')
    if record.get('source_hashes')!=source_hashes():
        raise ValueError('Current sources changed')
    for name,reviewer in REVIEWERS.items():
        compare_records(reviewer(),record['reviews'][name],name)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit('Output exists; choose a fresh path')
    record=generate();validate(record)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text());validate(frozen);compare_records(record,frozen)
    args.output.write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
    print(STATUS)
    print(f"{len(record['gates'])} gates pass; {len(record['inherited_provenance'])} inherited hashes agree.")
    print('Primary physical bounds only; nonlinear source classes, complete covariance means, full residual and blow-up remain open.')


if __name__=='__main__':
    main()
