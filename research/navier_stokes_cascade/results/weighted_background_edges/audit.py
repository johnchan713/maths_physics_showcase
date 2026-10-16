#!/usr/bin/env python3
"""Replay all-order weighted edges and the relative background stress estimate."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from edges import load,exact_review,derivative_review,diagonal_review,relative_edge_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='all-order-weighted-edges-and-relative-background-stress-audited'
FLAGS=dict(leading_all_fixed_order_weighted_edge_certificate_joined=True,
           every_relative_stress_coefficient_smooth_at_both_edges=True,
           additional_relative_and_weighted_diagonal_requirements_finite=True,
           actual_background_relative_stress_comparison_derived=True,
           actual_background_growth_and_oriented_cone_on_common_small_q_domain=True,
           relies_on_replayed_core_loop_repair_support_and_heat_lemmas=True,
           first_two_reserved_mean_pressure_tails_retained=True,
           actual_giant_profile_or_Borel_scales_numerically_evaluated=False,
           unconditional_uniform_physical_pulse_estimate_verified=False,
           complete_covariance_mean_update_verified=False,
           unconditional_complete_residual_improvement_verified=False,
           blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['verify_actual_physical_pulse_values_derivatives_and_supports_on_every_scale',
           'retain_first_two_pressure_tails_and_full_axial_velocity_in_complete_covariance_means',
           'prove_complete_residual_gain_with_viscosity_cutoffs_and_all_interactions']
DEPENDENCIES=('coefficient_support_audit','core_contraction_audit','core_growth_activation',
              'inner_collar_activation','outer_edge_stability','regional_cone_join',
              'background_jet_foundation','loop_modulation_bounds','repair_state_bounds')
REVIEWERS=dict(exact_flat_algebra_and_support_ledger=exact_review,
               both_edge_weights_and_fixed_derivatives=derivative_review,
               all_smaller_q_diagonal_schedule=diagonal_review,
               smooth_relative_heat_edge_and_direction=relative_edge_review)


def protocol():
    return dict(parent_commit='ede56e1149c4ea3728f1376039755a951f6fd71f',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                same_Md64_h_N_and_width=True,inner_exponent='t1^2/ya^2',outer_exponent='4/yb^2',
                exact_weight='exp(-t1^2/ya^2-4/yb^2)',distance='min(1,ya,yb)',
                derivative_family='Euler q*d_q at fixed X,eta; partial_X; partial_eta',
                total_order_weight_loss='3*m+3',
                additional_free_cutoff_requirements=['Tn/zeta jets through n','Tn/abs(T0) jets through n'],
                actual_cutoff_constants_or_scales_numerically_evaluated=False,
                relative_estimate='D^I[(That-T0)/abs(T0)]=O_I(q^(2h)) on the closed active annulus',
                normalized_stress='That=q^(A+1/2)*Tphysical',
                common_q_domain_for_every_fixed_jet=True,
                cone_domain_uses_only_value_and_first_velocity_jet_comparison=True,
                finite_N_not_reselected_for_higher_derivatives=True,
                mean_pressure_orders_one_and_two_not_discarded=True,
                diagnostic_h='1/200',diagnostics_are_not_continuum_proofs=True)


@lru_cache(None)
def parent_module():
    p=HERE.parent/'coefficient_support_audit'
    sys.path.insert(0,str(p))
    try:
        return load('weighted_parent_support_audit',p/'audit.py')
    finally:
        sys.path.pop(0)


def parent_review():
    module=parent_module();record=module.generate();module.validate(record)
    frozen=json.loads((HERE.parent/'coefficient_support_audit/evidence.json').read_text())
    module.compare_records(record,frozen)
    inherited=dict(record['inherited_provenance'])
    for path,expected in frozen['source_hashes'].items():
        inherited['coefficient_support_audit:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==expected
    return record,inherited


def source_hashes():
    paths=[HERE/n for n in ('README.md','edges.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():
        raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review()
    reviews={k:f() for k,f in REVIEWERS.items()}
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews.items()}
    gates['exact_parent_supports_replayed']=all(parent['gates'].values())
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
    if set(record.get('reviews',{}))!=set(REVIEWERS) or set(record.get('gates',{}))!=set(REVIEWERS)|{'exact_parent_supports_replayed','frozen_provenance_preserved'}:
        raise ValueError('Review or gate removed')
    if not all(record['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in record['reviews'].values()):
        raise ValueError('Weighted edge or relative comparison check failed')
    if len(record.get('inherited_provenance',{}))!=564 or not all(record['inherited_provenance'].values()):
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
    print('Analytic lemma dependencies retained; physical pulses, full covariance means, full residual and blow-up remain open.')


if __name__=='__main__':
    main()
