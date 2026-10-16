#!/usr/bin/env python3
"""Current realized field: conditional compact smooth force and preterminal energy."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from forcing import (load,domain_review,localization_review,terminal_review,
    extension_review,energy_review,admissibility_review)
from current_checks import diagnostic_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='conditional-current-field-localization-smooth-compact-force-and-energy'
FLAGS=dict(current_realized_complete_field_and_punctured_endpoint_jets_are_input=True,
    conditional_current_axisymmetric_localization_and_zero_initial_velocity_supplied=True,
    all_twelve_curl_pressure_viscous_and_nonlinear_terms_retained=True,
    force_divergence_in_exact_pressure_equation_retained=True,
    complete_current_terminal_force_jets_smooth_compact_and_flat_at_origin_conditionally_derived=True,
    conditional_one_smooth_compact_global_force_matching_exact_preterminal_residual_supplied=True,
    conditional_official_smooth_data_and_rapid_force_decay_requirements_met=True,
    conditional_current_whole_space_energy_and_integrated_dissipation_bound_supplied=True,
    conditional_complete_core_lower_bound_survives_current_localization=True,
    same_admitted_charts_parameters_width_and_Gaussian_exclusion_preserved=True,
    inherited_leading_background_inverse_and_uniform_physical_estimates_remain_analytical_premises=True,
    historical_frozen_records_reproduced_without_replacing_current_inputs=True,
    diagnostics_not_actual_global_fields_or_continuum_certificates=True,
    independent_full_analytic_foundation_review_complete=False,
    actual_terminal_force_norms_or_extension_scales_numerically_evaluated=False,
    raw_force_Taylor_convergence_or_uniform_terminal_order_constants_assumed=False,
    smooth_velocity_or_pressure_after_t_one_supplied=False,
    force_replaced_by_unforced_divergence_free_projection=False,
    Gaussian_cutoff_tails_forward_inverted=False,
    arbitrary_global_smooth_finite_energy_competitor_exclusion_verified=False,
    unconditional_admissible_force_and_solution_verified=False,
    blowup_verified=False,formal_proof_assistant_certificate=False)
REMAINING=['join_pressure_recovery_and_expanding_ball_comparison_to_same_current_localized_field_force_and_zero_datum',
    'control_all_spatial_infinity_fluxes_for_arbitrary_smooth_finite_energy_competitor_without_assumed_pressure_decay',
    'independently_audit_leading_profile_background_all_order_inverse_and_uniform_physical_bounds_before_unconditional_claims']
DEPENDENCIES=('current_endpoint_jets','current_potential_realization','finite_stage_induction',
    'localized_force_extension','full_potential_realization','physical_pulse_transfer',
    'compact_mean_corrections','background_jet_foundation')


def protocol():
    return dict(parent_commit='e4e38b2309d3a9d426c820b659df2a1bf80e15fc',
        paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        source_under_audit=['9.15-9.21','10.1-10.14'],
        official_definition='https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf; equations (4),(5),(7), alternative (C); checked 2026-10-16',
        actual_viscosity='1; other positive viscosities in manufactured identities only',
        current_input='complete current A,B,p0 and all fixed Cartesian residual/endpoint jets, not older hypothetical hierarchy',
        observation_domain='q_e=min(qstar/16,2^(-ell0-2)); tau0=q_e/16; z0=(q_e/16)^D; r0<R0',
        localized_field='u=curl(c*A)+c*B=c*u0+d; p=c*p0; d=grad(c) cross A; c axisymmetric',
        pressure='same full p_B tails and every wave/mean pressure difference; forced Poisson includes -div f',
        force_extension='sum_(j>=0) chi0(b_j*sigma)*sigma^j/j!*F_j; b0=2; increasing dyadic bj',
        extension_diagonal='spatial_order+time_order<=floor(j/2); derivative tail bounded by sum 2^-j',
        energy_bound='||u(t)||_2^2+2*integral_0^t ||grad u||_2^2<=F(t)^2; F(t)=integral_0^t ||f||_2',
        same_parameters_width_and_labels=True,
        diagnostic_field='manufactured complete current-endpoint potential/direct swirl/pressure on actual implicit q, including t=1 away from origin',
        diagnostic_positive_viscosities=[.6,1.,1.7],diagnostic_cutoff_scales=dict(radial=.55,axial=.4,time=.0015),
        actual_terminal_norms_scales_and_global_field_not_evaluated=True,
        scope='conditional current compact smooth force, zero datum, preterminal energy/dissipation and preserved core growth; not global competitor exclusion or unconditional blow-up')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'current_endpoint_jets';sys.path.insert(0,str(p))
    try:return load('current_force_endpoint_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    p=parent_module();r=p.generate();p.validate(r)
    frozen=json.loads((HERE.parent/'current_endpoint_jets/evidence.json').read_text())
    p.validate(frozen);p.compare_records(r,frozen)
    inherited=dict(r['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['current_endpoint_jets:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return r,inherited


@lru_cache(None)
def historical_review():
    # The current endpoint audit already reproduces this record. Reuse its
    # loaded module and explicitly review all old force/energy controls;
    # do not import the old field's conclusion as a current-field premise.
    parent_module().historical_review()
    p=sys.modules['current_endpoint_historical_localization']
    r=p.generate();p.validate(r)
    frozen=json.loads((HERE.parent/'localized_force_extension/evidence.json').read_text())
    p.validate(frozen);p.compare_records(r,frozen)
    controls={k:dict(checks={n:bool(v) for n,v in x['checks'].items()}) for k,x in r['reviews'].items()}
    return dict(old_force_field_not_substituted_for_current=True,controls=controls,
        checks=dict(all_historical_force_energy_and_terminal_controls_reproduce=all(r['gates'].values()),
            all_historical_source_bytes_unchanged=r['source_hashes']==frozen['source_hashes']))


REVIEWERS=dict(current_domain_and_support=domain_review,complete_localization_and_pressure=localization_review,
    global_terminal_force_jets=terminal_review,one_smooth_force_extension=extension_review,
    current_whole_space_energy_and_core_growth=energy_review,official_admissibility_interface=admissibility_review,
    current_q_localization_and_forced_pressure_controls=diagnostic_review,historical_force_controls_replayed=historical_review)
EXTRA_GATES={'current_actual_endpoint_and_realization_replayed','frozen_provenance_preserved'}


def source_hashes():
    paths=[HERE/n for n in ('README.md','forcing.py','current_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Frozen protocol changed')
    parent,inherited=parent_review();reviews={k:f() for k,f in REVIEWERS.items()}
    gates={k:bool(r['checks']) and all(r['checks'].values()) for k,r in reviews.items()}
    gates['current_actual_endpoint_and_realization_replayed']=all(parent['gates'].values()) and all(parent['flags'][k] for k in
        ('conditional_complete_potential_direct_swirl_and_pressure_terminal_jets_derived',
         'actual_flat_velocities_base_pressure_tails_and_Gaussian_exclusion_preserved',
         'current_realized_field_not_old_hypothetical_hierarchy_is_input'))
    gates['frozen_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=reviews,gates=gates,
        inherited_provenance=inherited,source_hashes=source_hashes(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():
        raise ValueError('Conditional current force, pressure, energy or scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Outstanding obligation changed')
    if set(r.get('reviews',{}))!=set(REVIEWERS) or set(r.get('gates',{}))!=set(REVIEWERS)|EXTRA_GATES:
        raise ValueError('Missing review or gate')
    if not all(r['gates'].values()) or not all(x.get('checks') and all(x['checks'].values()) for x in r['reviews'].values()):
        raise ValueError('Current force transfer gate failed')
    if len(r.get('inherited_provenance',{}))!=732 or not all(r['inherited_provenance'].values()):
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
    print('Conditional current compact force and energy only; competitor comparison, independent foundation and blow-up remain open.')


if __name__=='__main__':main()
