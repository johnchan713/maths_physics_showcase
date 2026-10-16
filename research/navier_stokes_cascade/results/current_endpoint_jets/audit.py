#!/usr/bin/env python3
"""Conditional current-field punctured terminal jets and a nonanalytic heat limit."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from endpoint import load,atlas_review,endpoint_review,heat_review,nonanalytic_review
from endpoint_checks import diagnostic_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='conditional-current-field-closed-time-atlas-and-punctured-terminal-jets'
FLAGS=dict(current_realized_field_not_old_hypothetical_hierarchy_is_input=True,
    closed_time_finite_band_and_mesh_geometry_audited=True,
    same_frozen_labels_absolute_pullbacks_and_physical_mean_bumps_retained=True,
    conditional_complete_potential_direct_swirl_and_pressure_terminal_jets_derived=True,
    conditional_uniform_Cauchy_and_finite_Taylor_remainder_bounds_derived=True,
    complete_exterior_heat_and_pressure_time_jets_derived=True,
    zero_radius_time_Taylor_series_for_nonzero_heat_exterior_and_pressure_proved=True,
    actual_flat_velocities_base_pressure_tails_and_Gaussian_exclusion_preserved=True,
    historical_frozen_evidence_reproduced_without_promotion=True,
    inherited_all_order_physical_background_and_inverse_bounds_remain_analytical_premises=True,
    diagnostics_are_not_actual_global_fields_or_continuum_certificates=True,
    independent_full_foundation_review_complete=False,
    actual_terminal_field_jets_or_infinite_coefficients_numerically_evaluated=False,
    fractional_auxiliary_map_proved_smooth_at_axis=False,
    analytic_or_zero_future_continuation_of_nonzero_terminal_heat_jets_justified=False,
    regular_velocity_or_pressure_at_singular_origin_inferred=False,
    current_complete_localization_initialization_admissible_force_and_energy_verified=False,
    unconditional_complete_residual_construction_verified=False,
    blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['join_complete_axisymmetric_physical_localization_and_zero_initialization_to_current_realized_field',
    'retain_all_curl_pressure_viscous_and_nonlinear_terms_and_extend_current_force_by_smooth_not_analytic_terminal_jet_gluing',
    'verify_whole_space_energy_dissipation_and_global_solution_comparison_for_same_current_field_and_force',
    'independently_review_leading_background_inverse_and_uniform_physical_bounds_before_unconditional_claims']
DEPENDENCIES=('current_potential_realization','finite_stage_induction','auxiliary_geometry_audit',
    'geometric_pulse_replay','physical_pulse_transfer','background_jet_foundation',
    'heat_exterior_audit','localized_force_extension')


def protocol():
    return dict(parent_commit='aac13c80650aa0f46320e723270431b9d2d64b3f',
        paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_under_audit=['6.1-6.22','9.15-9.21','10.1-10.14','Appendix A.7'],
        domain='same local cylinder; t<=1 away from x=0; no future field construction',
        endpoint_observation_domain='q_e=min(qstar/16,2^(-ell0-2)); tau0=q_e/16; z0=(q_e/16)^D; all active 4Q<qstar',
        representatives='complete A, separately divergence-safe direct swirl B, full seed pressure plus every correction difference',
        closed_band_window='ceil(log2(1/(2*qmax)))<=ell<=floor(log2(2/c)); ell>=ell0',
        closed_mesh='ceil(lo/delta)-1<=a<=floor(hi/delta)+1; T=0 uses a_T=-1,0,1',
        slow_chart_time='T=tau/Q; normalized time derivative coefficient -epsilon=-Q^h',
        current_cutoff_stage_bound='j<log2(qstar/c), since a_j>2^j/qstar',
        endpoint_input='uniform physical all-fixed-order bounds on q>=c, not interior smoothness or origin flatness alone',
        time_convention='tau=1-t; F_j=lim partial_t^j G; G(1-tau)=sum (-tau)^j F_j/j! plus finite remainder',
        heat_A='A_j=(h)_j*(1+h)_j',heat_S='S_j=sum binomial(j,k)*A_k*A_(j-k)',
        nonanalytic_condition='fixed h>0 and nonzero exterior constant C; smooth one-sided jets but zero Taylor radius',
        diagnostic_h_values=['1/200','9/1000'],diagnostic_q0='2/25',
        diagnostic_charts=['1/8','1/16'],diagnostic_tau=[0.0,0.000244140625,0.0009765625],
        diagnostic_jet_order=3,diagnostic_fields='manufactured potentials/direct swirl/pressure on actual implicit q; independent positive Gamma integrals',
        historical_records_preserved=True,parameters_width_geometry_and_primary_roots_unchanged=True,
        scope='conditional compatible punctured terminal jets; exact heat/pressure nonanalyticity; no current global-force or blow-up result')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'current_potential_realization';sys.path.insert(0,str(p))
    try:return load('current_endpoint_realization_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    p=parent_module();r=p.generate();p.validate(r)
    frozen=json.loads((HERE.parent/'current_potential_realization/evidence.json').read_text())
    p.validate(frozen);p.compare_records(r,frozen)
    inherited=dict(r['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['current_potential_realization:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return r,inherited


@lru_cache(None)
def geometry_review():
    g=load('current_endpoint_exact_geometry',HERE.parent/'auxiliary_geometry_audit/geometry.py')
    functions=dict(covering_gap=g.covering_gap_review,coloring=g.coloring_review,
        closed_separation=g.separation_review,common_copies=g.covering_review,counts=g.counting_review)
    return dict(historical_exact_geometry_not_a_pulse_amplitude_proof=True,
        checks={k:all(bool(v) for v in f()['checks'].values()) for k,f in functions.items()})


@lru_cache(None)
def historical_review():
    path=HERE.parent/'localized_force_extension';sys.path.insert(0,str(path))
    try:p=load('current_endpoint_historical_localization',path/'audit.py')
    finally:sys.path.pop(0)
    r=p.generate();p.validate(r);frozen=json.loads((path/'evidence.json').read_text())
    p.validate(frozen);p.compare_records(r,frozen)
    return dict(historical_force_and_energy_conclusions_not_imported_for_current_field=True,
        checks=dict(historical_record_and_all_diagnostics_reproduce=all(r['gates'].values()),
            inherited_endpoint_interface_uses_uniform_next_derivative=True,
            historical_sources_unchanged=r['source_hashes']==frozen['source_hashes']))


REVIEWERS=dict(finite_closed_time_atlas=atlas_review,complete_current_terminal_interface=endpoint_review,
    exact_exterior_swirl_and_pressure_jets=heat_review,nonanalytic_terminal_heat_limit=nonanalytic_review,
    manufactured_terminal_and_failure_controls=diagnostic_review,
    exact_global_geometry_replayed=geometry_review,historical_endpoint_transfer_replayed=historical_review)
EXTRA_GATES={'current_actual_realization_replayed','frozen_provenance_preserved'}


def source_hashes():
    paths=[HERE/n for n in ('README.md','endpoint.py','endpoint_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review();reviews={k:f() for k,f in REVIEWERS.items()}
    gates={k:bool(r['checks']) and all(r['checks'].values()) for k,r in reviews.items()}
    gates['current_actual_realization_replayed']=all(parent['gates'].values()) and all(parent['flags'][k] for k in
        ('conditional_single_real_smooth_divergence_free_local_field_for_q_positive_supplied',
         'conditional_complete_residual_flatness_in_every_fixed_Cartesian_jet_supplied',
         'designated_Gaussian_tails_never_inverted_or_independently_summed'))
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
        inherited_provenance=inherited,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():
        raise ValueError('Conditional endpoint scope, current seed or heat continuation changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Outstanding obligation changed')
    if set(r.get('reviews',{}))!=set(REVIEWERS) or set(r.get('gates',{}))!=set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(r['gates'].values()) or not all(v.get('checks') and all(v['checks'].values()) for v in r['reviews'].values()):
        raise ValueError('Endpoint gate failed')
    if len(r.get('inherited_provenance',{}))!=716 or not all(r['inherited_provenance'].values()):
        raise ValueError('Inherited provenance changed')
    if r.get('source_hashes')!=source_hashes():raise ValueError('Current sources changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')
    for k,f in REVIEWERS.items():compare_records(f(),r['reviews'][k],'root.reviews.'+k)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():raise SystemExit('Output exists; choose a fresh path')
    r=generate();validate(r)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text());validate(frozen);compare_records(r,frozen)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(STATUS)
    print(f"{len(r['gates'])} gates pass; {len(r['inherited_provenance'])} inherited hashes agree.")
    print('Conditional punctured terminal jets only; current global force, foundation review and blow-up remain open.')


if __name__=='__main__':main()
