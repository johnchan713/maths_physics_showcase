#!/usr/bin/env python3
"""Replay formal support identities without promoting physical pulse bounds."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from supports import load,geometry_review,exact_review,formal_review,inverse_review,source_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='formal-coefficient-supports-and-exact-reserved-mean-velocity-audited'
FLAGS=dict(actual_I3_I4_untouched_geometry_verified=True,
           finite_formal_coefficient_common_support_derived=True,
           exact_reserved_mean_velocity_identity_derived=True,
           first_two_pressure_tails_retained=True,
           outer_edge_common_support_premise_discharged=True,
           relies_on_analytic_collar_and_continuum_inverse_lemmas=True,
           all_order_weighted_leading_edge_certificate_verified=False,
           full_weighted_background_cutoff_inputs_verified=False,
           unconditional_uniform_actual_pulse_estimate_verified=False,
           complete_covariance_mean_update_verified=False,
           unconditional_complete_residual_improvement_verified=False,
           blowup_verified=False)
REMAINING=['audit_all_order_weighted_leading_and_first_stress_edges',
           'verify_all_scale_actual_pulse_values_derivatives_and_supports',
           'retain_first_two_pressure_tails_and_axial_velocity_in_complete_covariance_means',
           'prove_complete_residual_gain_with_viscosity_cutoffs_and_all_interactions']
DEPENDENCIES=('regional_cone_join','background_jet_foundation','background_summation',
              'annular_first_correction','outer_edge_stability','compact_mean_corrections')
REVIEWERS=dict(actual_reserved_geometry=geometry_review,exact_support_and_cutoff_identities=exact_review,
               finite_formal_induction=formal_review,continuum_linear_inverse=inverse_review,
               shifted_pressure_source_controls=source_review)


def protocol():
    return dict(parent_commit='3adeeaa46c64ee2dc6d47dcaf951baaf946cfc05',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                fixed_parameters_and_all_historical_evidence_unchanged=True,
                positive_velocity_edge='rho3*exp(3)',positive_bump_last_edge='rho3*exp(2.55)',
                reserved_mean_interval='rho3*exp((6,11))',
                exceptional_pressure_orders=[1,2],pressure_order3_and_above_edge='Xvel',
                raw_inverse_costs=['lambda^-1','e(eta)^-1'],
                full_lower_pressure_source_integrated_through_Xv=True,
                positive_order_velocity_cutoffs_are_potential_cutoffs=True,
                mean_pressure_not_equal_to_leading_pressure=True,
                actual_second_pressure_constant_not_evaluated=True,
                manufactured_h=.005,manufactured_profiles_not_actual_solution=True,
                finite_coefficients_not_uniform_in_order=True,
                inherited_analytic_lemmas_are_explicit_dependencies=True,
                no_full_residual_or_blowup_claim=True)


def parent_module():
    p=HERE.parent/'regional_cone_join'
    sys.path.insert(0,str(p))
    try:
        return load('support_parent_regional_audit',p/'audit.py')
    finally:
        sys.path.pop(0)


def parent_review():
    module=parent_module();record=module.generate();module.validate(record)
    frozen=json.loads((HERE.parent/'regional_cone_join/evidence.json').read_text())
    module.compare_records(record,frozen)
    inherited=dict(record['inherited_provenance'])
    for path,expected in frozen['source_hashes'].items():
        inherited['regional_cone_join:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==expected
    return record,inherited


def source_hashes():
    paths=[HERE/n for n in ('README.md','supports.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():
        raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review()
    reviews={k:f() for k,f in REVIEWERS.items()}
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews.items()}
    gates['parent_leading_join_replayed']=all(parent['gates'].values())
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
    if set(record.get('reviews',{}))!=set(REVIEWERS) or set(record.get('gates',{}))!=set(REVIEWERS)|{'parent_leading_join_replayed','frozen_provenance_preserved'}:
        raise ValueError('Review or gate removed')
    if not all(record['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in record['reviews'].values()):
        raise ValueError('Support check failed')
    if len(record.get('inherited_provenance',{}))!=551 or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited provenance changed')
    if record.get('source_hashes')!=source_hashes():
        raise ValueError('Current sources changed')
    # Scope and exceptional support rows must not rely on mutable success flags.
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
    print('Formal supports and exact mean velocity only; pressure tails, weighted edges and full residual limits retained.')


if __name__=='__main__':
    main()
