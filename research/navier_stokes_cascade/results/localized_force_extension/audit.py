#!/usr/bin/env python3
"""Conditional endpoint, full localization, smooth compact force and energy."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from endpoint_bounds import (borel_scale_review, domain_review, energy_ledger,
                             endpoint_interface_review)
from localization_checks import localization_review, support_axis_review
from force_checks import energy_review, heat_endpoint_review, terminal_gluing_review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'conditional-physical-localization-compact-force-and-whole-space-energy'
FLAGS = dict(conditional_endpoint_limits_from_inherited_uniform_physical_jets_supplied=True,
             conditional_complete_physical_localization_and_zero_initial_data_supplied=True,
             complete_curl_pressure_viscosity_and_transport_localization_terms_retained=True,
             conditional_all_order_compact_smooth_force_extension_supplied=True,
             conditional_whole_space_energy_and_integrated_dissipation_supplied=True,
             conditional_complete_localized_field_core_lower_bound_preserved=True,
             actual_global_profile_and_field_numerically_evaluated=False,
             actual_terminal_jets_or_extension_scales_evaluated=False,
             inherited_foundation_and_global_auxiliary_geometry_independently_verified=False,
             unconditional_admissible_force_and_solution_verified=False,
             arbitrary_global_smooth_bounded_energy_solution_comparison_verified=False,
             blowup_verified=False,
             diagnostic_success_is_not_a_continuum_proof=True,
             formal_proof_assistant_certificate=False)
REVIEW_NAMES = ('uniform_endpoint_interface', 'compact_localization_domain',
                'complete_localization_identity', 'axis_support_and_initialization',
                'heat_exterior_endpoint', 'all_order_force_extension_scales',
                'terminal_mixed_jet_gluing', 'whole_space_energy_and_dissipation')
GATE_NAMES = REVIEW_NAMES+('all_inherited_bytes_preserved', 'conditional_transfer_not_promoted')
REMAINING = ('arbitrary_global_smooth_bounded_energy_solution_comparison_with_pressure_and_spatial_infinity_fluxes',
             'independent_inherited_profile_background_all_order_physical_bounds_and_global_atlas_review',
             'only_after_these_inputs_global_solution_contradiction_and_unconditional_theorem_assessment')


def expected_protocol():
    return dict(status=STATUS,parent_commit='cef7a529973f9ff43d9107cd132663f6bfee7b6a',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['9.15-9.21','10.1-10.14'],
                actual_profile_parameters='unchanged inherited N, h, positive lambda, common qstar and global atlas hypotheses',
                actual_viscosity='1; other positive viscosities are algebra diagnostics only',
                localization_domain='tau0=qstar/16, z0=(qstar/16)^D; q<=3*qstar/8',
                terminal_extension='sum chi0(b_j*sigma)*sigma^j/j!*F_j; b_0=2, b_j increasing dyadic integers',
                extension_diagonal='spatial_order+time_order<=floor(j/2)',
                diagnostic_jet_points=[[1.12,-.4,-.19,.64],[1.3,.27,.25,.72],[1.44,-.9,.31,.8]],
                diagnostic_positive_viscosities=[.4,1.,1.7],diagnostic_wave_frequency=19,
                diagnostic_heat_h='1/200',diagnostic_heat_arguments=[0.,.3,2.],
                diagnostic_force_jet_terms=8,diagnostic_energy_quadrature=96,
                retained_failed_energy_quadrature=64,energy_integration_tolerance=2e-9,
                other_diagnostic_relative_tolerance=2e-10,
                diagnostic_data='manufactured full curls with Cartesian base, gamma integral, terminal jets, synthetic scale constants and smooth compact energy fixture',
                scope='Conditional endpoint limits, full physical localization, zero datum, smooth compact force, whole-space energy/dissipation and surviving core growth. Actual field evaluation, global-solution comparison and independent foundation/atlas review remain open.')


def proof_ledger():
    return dict(endpoint_uses_existing_uniform_physical_bounds_not_a_new_boundary_atlas=True,
                endpoint_Cauchy='|D_x^alpha partial_t^j G(t2)-D_x^alpha partial_t^j G(t1)|<=|t2-t1| sup |D_x^alpha partial_t^(j+1) G|',
                physical_endpoint_hypothesis_remains_conditional=True,
                heat_derivative_bound='|H^(j)(Z)|<=(h)_j*(1+h)_j for Z>=0',
                pressure_density_majorant='rho^(-3-4*h-2*j), integrable for rho>=r_min>0',
                cutoffs='c=chi_x(r^2,z)*chi_t(t); supported inside q<qstar/2',
                velocity='u=curl(c*A)+c*B=c*u0+d, d=grad c cross A',
                direct_swirl_divergence='partial_theta(c*b)=0; axisymmetry is required',
                full_localization_terms=['c*R0','c_t*u0','d_t','-nu*Delta(c)*u0','-2*nu*grad(c).grad(u0)',
                                         '-nu*Delta(d)','p0*grad(c)','c*(c-1)*u0.grad(u0)',
                                         'c*(u0.grad(c))*u0','c*u0.grad(d)','c*d.grad(u0)','d.grad(d)'],
                excluded_transport_is_identically_zero='d.grad(c)=0',
                pressure_Poisson='-Delta p=sum partial_i partial_j(u_i*u_j)-div f',
                force_jet_requirements='A:m+3, B:m+2, p0:m+1; one more derivative for uniform Cauchy',
                terminal_force_jets='F_j in C_c^infinity, support in K, every spatial jet at x=0 zero',
                extension_term_bound='C_jm*b_j^(m-j)*||F_j||_C^a',
                extension_diagonal='a+m<=floor(j/2); eventually includes each fixed mixed derivative',
                raw_terminal_Taylor_series_convergence_not_used=True,
                smooth_compact_force_support='K x [1-tau0,3/2]',
                no_regular_velocity_extension_through_t_one_inferred=True,
                energy='||u(t)||_2^2+2*integral_0^t ||grad u||_2^2<=F(t)^2, F(t)=integral_0^t ||f||_2',
                integrated_dissipation_finite_by_monotone_convergence=True,
                preserved_lower_bound='|u_theta(x_tau,1-tau)|>=(e0/2)*tau^(-A), for sufficiently small tau',
                arbitrary_global_solution_comparison_not_assumed=True)


def inherited_provenance():
    folders = ('paper_profile_audit','outer_pressure_pilot','axial_stress_audit',
               'intermediate_decay_audit','pulse_moment_audit','pulse_stress_audit',
               'post_pulse_stress_audit','heat_exterior_audit','axis_matching_audit',
               'axis_core_attachment','stress_realization_audit','joined_stress_construction',
               'angular_c2_transfer','compact_jet_envelope','loop_modulation_bounds',
               'repair_state_bounds','physical_residual_budget','first_correction_construction',
               'higher_jet_budget','annular_first_correction','background_summation',
               'physical_wave_interface','actual_band_pulse_bounds','compact_mean_corrections',
               'residual_cycle_closure','full_potential_realization')
    checks={}
    for folder in folders:
        record=json.loads((HERE.parent/folder/'evidence.json').read_text())
        if 'source_sha256' in record:
            hashes,root=record['source_sha256'],HERE.parent/folder
        elif 'source_hashes' in record:
            hashes,root=record['source_hashes'],PROJECT
        else:
            hashes,root=record.get('provenance',{}),None
        if not hashes:
            raise ValueError('Missing inherited ledger: '+folder)
        for path,digest in hashes.items():
            routed=root if root is not None else (PROJECT if path.startswith('results/') or path=='GOAL.md' else HERE.parent)
            checks[folder+':'+path]=hashlib.sha256((routed/path).read_bytes()).hexdigest()==digest
    return checks


def source_hashes():
    names=('README.md','.gitignore','requirements.txt','protocol.json','endpoint_bounds.py',
           'localization_checks.py','force_checks.py','audit.py','test_audit.py')
    paths=[HERE/name for name in names]+[PROJECT/'GOAL.md',
          HERE.parent/'full_potential_realization/evidence.json',
          HERE.parent/'heat_exterior_audit/evidence.json']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def encode(value):
    if isinstance(value,Fraction):
        return str(value)
    if isinstance(value,dict):
        return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [encode(v) for v in value]
    if isinstance(value,np.bool_):
        return bool(value)
    if isinstance(value,(float,np.floating)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite diagnostic')
        return float(value)
    return value


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != expected_protocol():
        raise ValueError('Frozen scope or parameters changed')
    energy=energy_review()
    energy['analytic_ledger']=energy_ledger()
    energy['checks'].update(energy['analytic_ledger']['checks'])
    reviews=dict(uniform_endpoint_interface=endpoint_interface_review(),
                 compact_localization_domain=domain_review(),
                 complete_localization_identity=localization_review(),
                 axis_support_and_initialization=support_axis_review(),
                 heat_exterior_endpoint=heat_endpoint_review(),
                 all_order_force_extension_scales=borel_scale_review(),
                 terminal_mixed_jet_gluing=terminal_gluing_review(),
                 whole_space_energy_and_dissipation=energy)
    inherited=inherited_provenance()
    gates={name:all(reviews[name]['checks'].values()) for name in REVIEW_NAMES}
    gates['all_inherited_bytes_preserved']=bool(inherited) and all(inherited.values())
    gates['conditional_transfer_not_promoted']=not any(FLAGS[k] for k in
        ('actual_global_profile_and_field_numerically_evaluated','actual_terminal_jets_or_extension_scales_evaluated',
         'inherited_foundation_and_global_auxiliary_geometry_independently_verified',
         'unconditional_admissible_force_and_solution_verified',
         'arbitrary_global_smooth_bounded_energy_solution_comparison_verified','blowup_verified'))
    return encode(dict(status=STATUS,protocol=expected_protocol(),flags=FLAGS,
                       conditional_proof=proof_ledger(),reviews=reviews,
                       inherited_provenance=inherited,source_hashes=source_hashes(),
                       gates=gates,remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS or record.get('conditional_proof') != proof_ledger():
        raise ValueError('Conditional forcing transfer is not an unconditional blow-up construction')
    if record.get('protocol') != expected_protocol() or record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('Frozen scope or remaining obligations changed')
    if set(record.get('gates',{})) != set(GATE_NAMES) or not all(record['gates'].values()):
        raise ValueError('A required localization, endpoint, extension, energy or provenance gate failed')


def compare_records(actual,frozen,location='root'):
    if isinstance(actual,dict):
        if not isinstance(frozen,dict) or set(actual) != set(frozen):
            raise ValueError('Frozen fields differ at '+location)
        for key in actual:
            compare_records(actual[key],frozen[key],location+'.'+key)
    elif isinstance(actual,list):
        if not isinstance(frozen,list) or len(actual) != len(frozen):
            raise ValueError('Frozen list differs at '+location)
        for i,(a,b) in enumerate(zip(actual,frozen)):
            compare_records(a,b,location+'['+str(i)+']')
    elif location.startswith('root.reviews.') and isinstance(actual,float):
        if isinstance(frozen,bool) or not isinstance(frozen,(float,int)) or not math.isclose(actual,frozen,rel_tol=1e-7,abs_tol=2e-10):
            raise ValueError('Frozen diagnostic differs at '+location)
    elif type(actual) is not type(frozen) or actual != frozen:
        raise ValueError('Frozen exact arithmetic or metadata differs at '+location)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit('Output already exists; choose a fresh path')
    record=generate()
    validate(record)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text())
        validate(frozen)
        compare_records(record,frozen)
    args.output.write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
    print(STATUS)
    print(str(len(record['gates']))+' gates pass; endpoint, compact force and whole-space energy remain conditional.')
    print('Global-solution comparison and independent foundation/atlas review remain open.')


if __name__=='__main__':
    main()
