#!/usr/bin/env python3
"""Constructive auxiliary support geometry; no unconditional blow-up claim."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from geometry import (coloring_review,counting_review,covering_gap_review,
                       covering_review,separation_review)
from checks import (chain_rule_review,compatibility_review,failed_shortcuts_review,
                     partition_review,path_review,scale_and_source_review)

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='constructive-all-band-auxiliary-supports-and-common-chart-geometry'
FLAGS=dict(explicit_all_band_coloring_centers_and_closed_separation_certificate_supplied=True,
           exact_covering_counts_and_Haar_compatibility_supplied=True,
           geometric_pullback_and_liftwise_path_compatibility_supplied=True,
           source_fast_scales_and_actual_physical_chain_rule_rederived=True,
           source_outer_parentheses_confirmed_by_visual_PDF_inspection=True,
           suspected_factor_two_normalization_error_rejected=True,
           prior_frozen_diagnostic_formula_matches_the_rendered_source=True,
           full_global_pulse_and_nonlinear_correction_hierarchy_revalidated=False,
           leading_profile_background_and_all_order_inverse_foundation_independently_verified=False,
           actual_global_profile_and_field_numerically_evaluated=False,
           unconditional_blowup_construction_verified=False,blowup_verified=False,
           diagnostic_success_is_not_a_continuum_proof=True,
           formal_proof_assistant_certificate=False)
REVIEW_NAMES=('rendered_source_and_floor_scales','all_band_covering_gap','slow_partition_coverage',
              'explicit_global_coloring','rational_centers_and_closed_separation','common_covering_copies_and_Haar',
              'pointwise_and_radial_label_counts','lifted_pulse_paths','actual_physical_chain_rule',
              'global_pullback_compatibility','failed_geometric_shortcuts')
GATE_NAMES=REVIEW_NAMES+('all_inherited_bytes_preserved','geometry_not_promoted_to_full_construction')
REMAINING=('primary_pulse_phase_frame_and_source_estimates_with_the_explicit_fixed_geometry',
           'direct_leading_profile_background_collar_mean_inverse_and_parameter_audit',
           'full_global_nonlinear_assembly_review_before_unconditional_exclusion_assessment')


def expected_protocol():
    return dict(status=STATUS,parent_commit='464244c7334eb0a73579431e5155f1699565f61a',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['6.1-6.22','Lemmas 6.1-6.3'],
                visual_source_check='PDF page 63 equation (6.2); outer 2 multiplies both terms',
                canonical_radial_exponent='2*((1+h)*rho_g-h*kappa_s)',
                rejected_text_extraction_variant='2*(1+h)*rho_g-h*kappa_s',
                minimum_band=16,neighboring_band_difference=4,maximum_covering_difference=3,
                tensor_mesh='delta_ell=ell^-6=S_ell^-3; original halfwidth delta, doubled halfwidth 2*delta',
                colors=1250,color_period=5,prime_denominator=2003,
                centers='c_n=(0,n/2003), n=1,...,1250',radius_upper_bound='1/6057072',
                center_congruence_margin='3/8012; not a uniform absolute-torus Euclidean distance',
                maximum_common_torus_copies=2744,maximum_common_transition_matrix_norm=188,
                original_pointwise_slow_overlap_bound=162,
                diagnostic_h_values=['1/200','9/1000'],
                diagnostic_bands=[16,17,32,160,1600,1000000,1000000000,1000000000000],
                diagnostic_jet_bands=[16,20,24],diagnostic_relative_tolerance=2e-12,
                diagnostic_data='exact integer/rational geometry, high precision floor/path identities and manufactured scalar evaluated on the actual absolute map',
                scope='Constructive pure support/covering geometry and preservation of compatibility by the specified physical operators. Full pulse/background/nonlinear estimates and unconditional blow-up remain unverified.')


def proof_ledger():
    return dict(source_formula='d_r=2*((1+h)*rho-h*kappa); visual check rejects a missing-parentheses false alarm',
                exact_fast_scales='c_i=T_g^(-theta)*S^-1; M_i=Lambda_g^(-theta)*epsilon^(-kappa)*S^(-rho), theta in [0,1)',
                all_band_covering_gap='|i(ell)-i(ell_prime)|<=1873/625<3 when ell,ell_prime>=16 and |ell-ell_prime|<=4',
                global_color='1+sign+2*(a1%5+5*(a2%5+5*(a3%5+5*(ell%5))))',
                color_count=1250,centers='(0,n/2003)',
                covering_off_diagonal_entries=[1,8,50],
                exact_separation='center gap>=1/2003, enlarged displacement<=1/8012 for r0<=1/6057072',
                absolute_preimages_disjoint_without_claiming_uniform_absolute_distance=True,
                graph_and_partition_include_closed_supports_and_tau_zero_before_label_selection=True,
                signs_duplicate_cutoff_but_squared_partition_is_counted_once_per_box=True,
                lift_count='14^Delta<=2744',common_matrix_norm_bound=188,
                pointwise_slow_label_bound=162,compact_box_label_count='O(S^9)',radial_line_label_count='O(S^3)',
                radial_integration_leaves_q_and_band_set_fixed=True,
                absolute_phase_map='Y=v_r*r^d_r+v_t*t mod Z^2',
                physical_operators='Dr=partial_R+M_i*d_r*R^(d_r-1)*L_i; Dz=epsilon*partial_Z; tstar=-epsilon*partial_T+c_i*N_i',
                pulse_path='H_w=H+T_g^(-Delta)*c_i*(w-v(H))*v_t',
                pulse_path_derivative='I-v_t*lambda_t; higher derivatives zero',
                compatibility='define once on absolute Y; fixed physical radial/mean operators and liftwise ODE uniqueness preserve pullbacks',
                independent_fixed_R_chart_bumps_or_averaging_across_copies_not_allowed=True,
                smooth_axis_extension_comes_from_supported_fields_not_fractional_map_smoothness=True,
                analytic_leading_profile_and_global_PDE_hierarchy_not_reproved=True)


def inherited_provenance():
    folders=('paper_profile_audit','outer_pressure_pilot','axial_stress_audit','intermediate_decay_audit',
             'pulse_moment_audit','pulse_stress_audit','post_pulse_stress_audit','heat_exterior_audit',
             'axis_matching_audit','axis_core_attachment','stress_realization_audit','joined_stress_construction',
             'angular_c2_transfer','compact_jet_envelope','loop_modulation_bounds','repair_state_bounds',
             'physical_residual_budget','first_correction_construction','higher_jet_budget','annular_first_correction',
             'background_summation','physical_wave_interface','actual_band_pulse_bounds','compact_mean_corrections',
             'residual_cycle_closure','full_potential_realization','localized_force_extension','global_solution_comparison')
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
    names=('README.md','.gitignore','requirements.txt','protocol.json','geometry.py','checks.py','audit.py','test_audit.py')
    paths=[HERE/name for name in names]+[PROJECT/'GOAL.md',HERE.parent/'global_solution_comparison/evidence.json',
                                       HERE.parent/'actual_band_pulse_bounds/evidence.json']
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
    if json.loads((HERE/'protocol.json').read_text())!=expected_protocol():
        raise ValueError('Frozen scope or parameters changed')
    reviews=dict(rendered_source_and_floor_scales=scale_and_source_review(),
                 all_band_covering_gap=covering_gap_review(),slow_partition_coverage=partition_review(),
                 explicit_global_coloring=coloring_review(),rational_centers_and_closed_separation=separation_review(),
                 common_covering_copies_and_Haar=covering_review(),pointwise_and_radial_label_counts=counting_review(),
                 lifted_pulse_paths=path_review(),actual_physical_chain_rule=chain_rule_review(),
                 global_pullback_compatibility=compatibility_review(),failed_geometric_shortcuts=failed_shortcuts_review())
    inherited=inherited_provenance()
    gates={name:all(reviews[name]['checks'].values()) for name in REVIEW_NAMES}
    gates['all_inherited_bytes_preserved']=bool(inherited) and all(inherited.values())
    gates['geometry_not_promoted_to_full_construction']=not any(FLAGS[k] for k in
        ('full_global_pulse_and_nonlinear_correction_hierarchy_revalidated',
         'leading_profile_background_and_all_order_inverse_foundation_independently_verified',
         'actual_global_profile_and_field_numerically_evaluated','unconditional_blowup_construction_verified','blowup_verified'))
    return encode(dict(status=STATUS,protocol=expected_protocol(),flags=FLAGS,geometric_proof=proof_ledger(),
                       reviews=reviews,inherited_provenance=inherited,source_hashes=source_hashes(),
                       gates=gates,remaining_obligations=REMAINING))


def validate(record):
    if record.get('status')!=STATUS or record.get('flags')!=FLAGS or record.get('geometric_proof')!=proof_ledger():
        raise ValueError('Pure geometric construction does not verify the full analytic hierarchy or blow-up')
    if record.get('protocol')!=expected_protocol() or record.get('remaining_obligations')!=list(REMAINING):
        raise ValueError('Frozen scope or remaining obligations changed')
    if set(record.get('gates',{}))!=set(GATE_NAMES) or not all(record['gates'].values()):
        raise ValueError('A required geometry, compatibility, scope or provenance gate failed')


def compare_records(actual,frozen,location='root'):
    if isinstance(actual,dict):
        if not isinstance(frozen,dict) or set(actual)!=set(frozen):
            raise ValueError('Frozen fields differ at '+location)
        for key in actual:
            compare_records(actual[key],frozen[key],location+'.'+key)
    elif isinstance(actual,list):
        if not isinstance(frozen,list) or len(actual)!=len(frozen):
            raise ValueError('Frozen list differs at '+location)
        for i,(a,b) in enumerate(zip(actual,frozen)):
            compare_records(a,b,location+'['+str(i)+']')
    elif location.startswith('root.reviews.') and isinstance(actual,float):
        if isinstance(frozen,bool) or not isinstance(frozen,(float,int)) or not math.isclose(actual,frozen,rel_tol=1e-7,abs_tol=2e-10):
            raise ValueError('Frozen diagnostic differs at '+location)
    elif type(actual)is not type(frozen) or actual!=frozen:
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
    print(str(len(record['gates']))+' gates pass; full pulse/profile/PDE estimates and blow-up remain unverified.')
    print('Visual source check rejected, rather than confirmed, the suspected normalization error.')


if __name__=='__main__':
    main()
