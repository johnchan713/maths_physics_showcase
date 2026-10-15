#!/usr/bin/env python3
"""General finite-energy comparison and a conditional global-exclusion step."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from comparison_bounds import (flux_review,gronwall_review,holder_review,kernel_review,
                               missing_hypothesis_review,pressure_recovery_review,
                               young_absorption_review)
from pressure_checks import (commutator_review,fourier_pressure_review,
                             negative_sobolev_review,weighted_energy_review)

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='general-finite-energy-comparison-and-conditional-global-solution-exclusion'
FLAGS=dict(general_smooth_compact_solution_comparison_proof_supplied=True,
           pressure_gradient_recovered_without_assumed_pressure_growth_or_decay=True,
           no_global_derivative_or_dissipation_bound_on_competitor_assumed=True,
           radius_independent_flux_constants_and_vanishing_remainder_supplied=True,
           conditional_global_smooth_bounded_energy_solution_exclusion_supplied=True,
           all_downstream_steps_remain_subject_to_inherited_continuum_inputs=True,
           actual_global_profile_and_field_numerically_evaluated=False,
           inherited_foundation_and_global_auxiliary_geometry_independently_verified=False,
           unconditional_admissible_construction_and_exclusion_verified=False,
           unforced_blowup_verified=False,blowup_verified=False,
           diagnostic_success_is_not_a_continuum_proof=True,
           periodic_diagnostics_are_not_whole_space_infinity_estimates=True,
           formal_proof_assistant_certificate=False)
REVIEW_NAMES=('pressure_recovery_and_negative_Sobolev','complete_difference_pressure_algebra',
              'commutator_kernel_and_independent_symbols','weighted_Holder_and_cutoff_powers',
              'complete_weighted_energy_identity','flux_absorption_and_radius_decay',
              'sharp_failed_shortcuts','fixed_T_comparison_and_conditional_exclusion')
GATE_NAMES=REVIEW_NAMES+('all_inherited_bytes_preserved','conditional_route_not_promoted')
REMAINING=('global_auxiliary_atlas_and_support_foundation_audit_from_defining_maps',
           'independent_leading_profile_background_all_order_inverse_and_parameter_review',
           'unconditional_theorem_assessment_only_after_all_inherited_continuum_inputs_are_verified')


def expected_protocol():
    return dict(status=STATUS,parent_commit='b510fd042c9f41709b5e685cd84ee865352e5f31',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['10.15-10.21','Lemma 10.5'],
                standard_external_input='classical Riesz transforms bounded on Lp, 1<p<infinity; applied twice at p=3/2',
                primary_standard_input_reference='https://www.math.stonybrook.edu/~aknapp/download/a2-1-realanal-clickable.pdf',
                primary_reference_location='Knapp, Advanced Real Analysis, Chapter III sections 5-6',
                actual_profile_parameters='unchanged; comparison lemma does not use pulse-specific parameters',
                actual_viscosity='1',competitor_hypotheses='smooth on every compact set; uniform L2 on each fixed [0,T]; same force and datum',
                no_pressure_growth_or_decay_assumption=True,no_global_competitor_gradient_assumption=True,
                cutoff_power=8,commutator_Lp='4/3',commutator_radius_power='-3/4',
                absorbed_energy_remainder_radius_power='-4/3',diagnostic_periodic_grid=32,
                diagnostic_intervals='one instantaneous PDE jet; no global flow trajectory or infinity solve',
                diagnostic_radial_radii=[1,3,19],diagnostic_energy_radii=[2,8,32,128],
                diagnostic_float_relative_tolerance=2e-12,
                scope='General finite-energy comparison lemma and conditional exclusion for the inherited localized field. Foundation/global-atlas validity remains open; periodic algebra checks do not certify whole-space bounds or unconditional blow-up.')


def proof_ledger():
    return dict(tensor='g=v tensor v-u tensor u=w tensor w+w tensor u+u tensor w in L_infinity_t L1_x',
                difference_equation='w_t+v.grad w+w.grad u=Delta w-grad pi; div w=0',
                pressure_multiplier='pi_star_hat=-sum xi_i*xi_j/|xi|^2*g_ij_hat',
                pressure_in_negative_Sobolev='pi_star in L_infinity_t H^-2; averaged pressure gradients in H^-3',
                harmonic_ambiguity_removed_by='time averaged equation plus harmonic H^-3 Fourier-support argument',
                no_pressure_growth_or_decay_and_no_global_w_time_derivative_norm_assumed=True,
                cutoff='phi_R^8, phi_R=phi(x/R), phi=1 on unit ball and zero outside radius two',
                weighted_Sobolev='B_R=||phi_R^4*w||_6<=C*(D_R+R^-1*W_T)',
                pressure_split='phi_R^4*pi_star=sum T_ij(phi_R^4*g_ij)+sum [phi_R^4,T_ij]g_ij',
                commutator_kernel='|x|^-3*min(|x|/R,1), L^(4/3) norm <=C*R^(-3/4)',
                flux='pressure<=C_T*(R^-1*(B_R^(3/2)+B_R^(1/2))+R^(-7/4)*B_R^(3/4))',
                transport_support='v=w on support grad chi_R once R exceeds support radius of u',
                young_remainder_powers=['4','4/3','14/5'],
                final_comparison_inequality='0.5*E_R_prime+0.5*D_R^2<=||grad u||_infinity*E_R+C_T*R^(-4/3)',
                every_constant_in_final_inequality_is_radius_independent=True,
                radius_limit_after_fixed_T='E_R(t)<=K_T*R^(-4/3), so w=0 on each fixed ball',
                conditional_exclusion='a globally smooth competitor equals u before 1, contradicting growth along points approaching (0,1)',
                independent_profile_or_atlas_validation_not_supplied=True,
                next_gate_is_foundation_audit_not_more_downstream_conditional_transfer=True)


def inherited_provenance():
    folders=('paper_profile_audit','outer_pressure_pilot','axial_stress_audit','intermediate_decay_audit',
             'pulse_moment_audit','pulse_stress_audit','post_pulse_stress_audit','heat_exterior_audit',
             'axis_matching_audit','axis_core_attachment','stress_realization_audit','joined_stress_construction',
             'angular_c2_transfer','compact_jet_envelope','loop_modulation_bounds','repair_state_bounds',
             'physical_residual_budget','first_correction_construction','higher_jet_budget','annular_first_correction',
             'background_summation','physical_wave_interface','actual_band_pulse_bounds','compact_mean_corrections',
             'residual_cycle_closure','full_potential_realization','localized_force_extension')
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
    names=('README.md','.gitignore','requirements.txt','protocol.json','comparison_bounds.py',
           'pressure_checks.py','audit.py','test_audit.py')
    paths=[HERE/name for name in names]+[PROJECT/'GOAL.md',HERE.parent/'localized_force_extension/evidence.json']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def combine(**reviews):
    checks={}
    for name,review in reviews.items():
        checks.update({name+':'+key:value for key,value in review['checks'].items()})
    return dict(reviews=reviews,checks=checks)


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
    reviews=dict(pressure_recovery_and_negative_Sobolev=combine(analytic=pressure_recovery_review(),radial=negative_sobolev_review()),
                 complete_difference_pressure_algebra=fourier_pressure_review(),
                 commutator_kernel_and_independent_symbols=combine(kernel=kernel_review(),symbols=commutator_review()),
                 weighted_Holder_and_cutoff_powers=holder_review(),
                 complete_weighted_energy_identity=weighted_energy_review(),
                 flux_absorption_and_radius_decay=combine(exponents=flux_review(),maximization=young_absorption_review()),
                 sharp_failed_shortcuts=missing_hypothesis_review(),
                 fixed_T_comparison_and_conditional_exclusion=gronwall_review())
    inherited=inherited_provenance()
    gates={name:all(reviews[name]['checks'].values()) for name in REVIEW_NAMES}
    gates['all_inherited_bytes_preserved']=bool(inherited) and all(inherited.values())
    gates['conditional_route_not_promoted']=not any(FLAGS[k] for k in
        ('actual_global_profile_and_field_numerically_evaluated',
         'inherited_foundation_and_global_auxiliary_geometry_independently_verified',
         'unconditional_admissible_construction_and_exclusion_verified','unforced_blowup_verified','blowup_verified'))
    return encode(dict(status=STATUS,protocol=expected_protocol(),flags=FLAGS,conditional_proof=proof_ledger(),
                       reviews=reviews,inherited_provenance=inherited,source_hashes=source_hashes(),
                       gates=gates,remaining_obligations=REMAINING))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS or record.get('conditional_proof') != proof_ledger():
        raise ValueError('Conditional comparison/exclusion does not verify the inherited foundation or unconditional blow-up')
    if record.get('protocol') != expected_protocol() or record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('Frozen scope or remaining obligations changed')
    if set(record.get('gates',{})) != set(GATE_NAMES) or not all(record['gates'].values()):
        raise ValueError('A required pressure, flux, comparison, scope or provenance gate failed')


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
    print(str(len(record['gates']))+' gates pass; the inherited construction and global exclusion remain conditional.')
    print('Next: audit the global auxiliary atlas/support and the inherited continuum foundation.')


if __name__=='__main__':
    main()
