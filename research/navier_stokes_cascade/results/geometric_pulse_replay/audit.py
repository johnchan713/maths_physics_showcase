#!/usr/bin/env python3
"""Conditional radius-explicit pulse replay; no unconditional field claim."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
import numpy as np
from estimates import (RADIUS, derivative_and_tail_review, exact_scale_review,
                       threshold_review)
from diagnostics import (covariance_review, failed_shortcuts_review, frame_review,
                         harmonic_factor_review, lifted_source_review, sourced_review)

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='conditional-pulse-and-source-replay-with-explicit-global-geometry'
FLAGS=dict(radius_explicit_phase_frame_and_positive_pulse_argument_supplied=True,
           transverse_area_and_actual_Haar_covariance_radius_cost_retained=True,
           common_domain_source_inverse_and_fixed_order_jet_argument_supplied=True,
           full_cutoff_tail_and_liftwise_source_support_conditions_retained=True,
           manufactured_fixture_uses_globally_admissible_radius=True,
           old_single_band_radius_not_promoted_to_global_geometry=True,
           numerical_diagnostics_are_not_continuum_estimates=True,
           actual_global_leading_profile_or_summed_background_numerically_evaluated=False,
           actual_profile_value_gap_and_jet_constants_interval_certified=False,
           leading_profile_background_and_mean_inverse_foundation_independently_verified=False,
           full_global_nonlinear_hierarchy_independently_revalidated=False,
           unconditional_complete_residual_improvement_verified=False,
           unconditional_blowup_construction_verified=False,blowup_verified=False,
           formal_proof_assistant_certificate=False)
REVIEW_NAMES=('exact_radius_and_scale_powers','all_later_domain_threshold',
              'global_radius_phase_and_frame','full_Haar_primary_covariance',
              'independent_sourced_projected_ODE','exact_higher_harmonic_factor',
              'liftwise_mixed_band_sources','fixed_order_jets_and_flat_tails','failed_shortcuts')
GATE_NAMES=REVIEW_NAMES+('refined_ODE_and_full_interval_quadrature_agree',
                         'all_inherited_bytes_preserved','conditional_scope_not_promoted')
REMAINING=('direct_leading_profile_strict_gaps_edge_bounds_and_common_collar_at_fixed_N_h',
           'summed_background_all_fixed_jet_estimates_and_all_order_mean_inverse_foundation',
           'independent_global_nonlinear_and_cutoff_assembly_review_before_unconditional_exclusion')


def expected_protocol():
    return dict(status=STATUS,parent_commit='1d7a71fd249cdb136f2f8658f5fa85afb740f5dc',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['6.10-6.22','7.2-7.29','Proposition 7.2','Lemmas 7.1,7.4,7.5'],
                global_radius=str(RADIUS),diagnostic_radii=[str(RADIUS),str(RADIUS/2)],
                diagnostic_bands=[51600,103200],diagnostic_h='1/200',actual_h='fixed inherited h in (0,1/100)',
                actual_leading_radial_frequency='unchanged inherited N=1+floor(H^32)',
                diagnostic_precision=100,diagnostic_frame_entry_bound_times_length=30,
                diagnostic_frame_bound_is_sampled_not_interval_certified=True,
                transverse_cutoff='beta(2*xi/r0), closed support [-r0/2,r0/2] strictly inside (-r0,r0)',
                ODE_solver='Radau after division by exact reference P; separate 2-frame and 3-component sourced equations',
                covariance_quadrature='full pulse in (v-L/2)/sqrt(L), with relative pulse and transverse cutoffs',
                relative_refinement_tolerance=3e-8,maximum_independent_sourced_ODE_error=2e-8,
                actual_profile_constants_or_actual_qstar_numerically_evaluated=False,
                scope='Conditional radius-explicit pulse/source estimates; actual strict profile gaps, summed background jets and all-order mean inverse foundation remain unverified independently. No unconditional full residual improvement or blow-up.')


def proof_ledger():
    return dict(length='2*r0*S<=L=2*r0/ci<2*r0*Tg*S<S',
                phase_smallness='S^2*(epsilon+epsilon^2+k^-1)<=1; |p-ptilde|<=1/k with nonzero carrier',
                frame='entrywise |E|<=CE/L; |d-dref|<=Cd/L, constants independent of radius on accepted long pulses',
                positive_pulse='Riccati interval CE/(mu*L)<=1/8 and |log(zplus/P)|<=Cd+CE*(1+CE/(mu*L))',
                forward_inverse='||Vm(v,w)||<=exp(2*CE+Cd)*P(v)/P(w), all nonzero integer m',
                source='zero data; amplitude value cost L, pressure extra epsilon^1/2; fixed-order derivative costs polynomial in S',
                full_Haar_weight='h=(1+b_g^2)/2 * integral chi_g^2 dxi * ci * integral psi^2*x^2 dv',
                covariance='h comparable to r0^2/sqrt(L), hence r0^(3/2)/sqrt(S)',
                inverse='norm H^-1 <= C*r0^(-3/2)*sqrt(S); determinant >= c*r0^3/S',
                primary_root='r0^(-3/4)*S^(1/4); jth derivatives add S^(2*j) and retained edge losses',
                complete_primary_coefficient='C_I*r0^(-3/4)*epsilon^(1/2)*S^(3*|I|+1/4)*sqrt(zeta)*delta^(-b_I)*P',
                domain='first band above inherited preliminary band, sqrt(Lmin/(2*r0)) and all-later exponential threshold; Lmin includes every finite normal, frame and cone value gap',
                mixed_sources='solve on each common-torus lift; never average different preimages',
                temporal_support='no extension before psi; original global source requires smooth endpoint zero extension',
                tails='P<=exp(-cG*L/50)<=exp(-(cG*r0/25)*ell^2) outside central cutoff plateau',
                radius_does_not_add_epsilon_power_loss=True,
                actual_profile_and_all_order_background_hypotheses_retained=True,
                continuum_arguments_are_not_inferred_from_finite_diagnostics=True)


def inherited_provenance():
    folders=('paper_profile_audit','outer_pressure_pilot','axial_stress_audit','intermediate_decay_audit',
             'pulse_moment_audit','pulse_stress_audit','post_pulse_stress_audit','heat_exterior_audit',
             'axis_matching_audit','axis_core_attachment','stress_realization_audit','joined_stress_construction',
             'angular_c2_transfer','compact_jet_envelope','loop_modulation_bounds','repair_state_bounds',
             'physical_residual_budget','first_correction_construction','higher_jet_budget','annular_first_correction',
             'background_summation','physical_wave_interface','actual_band_pulse_bounds','compact_mean_corrections',
             'residual_cycle_closure','full_potential_realization','localized_force_extension','global_solution_comparison',
             'auxiliary_geometry_audit')
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
    names=('README.md','.gitignore','requirements.txt','protocol.json','estimates.py','diagnostics.py','audit.py','test_audit.py')
    paths=[HERE/name for name in names]+[PROJECT/'GOAL.md',HERE.parent/'auxiliary_geometry_audit/evidence.json',
                                       HERE.parent/'actual_band_pulse_bounds/evidence.json']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def encode(value):
    if isinstance(value,Fraction):
        return str(value)
    if isinstance(value,mp.mpf):
        if not mp.isfinite(value):
            raise ValueError('Nonfinite diagnostic')
        return float(value)
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


def refinement_review(coarse,refined):
    relative=0.
    for a,b in zip(coarse['rows'],refined['rows']):
        for key in ('h','h_times_sqrt_length_over_radius_squared','minimum_normalized_radial_pulse'):
            relative=max(relative,abs(a[key]-b[key])/max(abs(b[key]),1e-300))
        relative=max(relative,float(np.linalg.norm(np.array(a['column'])-b['column']))/np.linalg.norm(b['column']))
    return dict(maximum_relative_covariance_and_radial_refinement_error=relative,
                checks=dict(refinement_with_tighter_ODE_and_full_interval_quadrature=relative<expected_protocol()['relative_refinement_tolerance']))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=expected_protocol():
        raise ValueError('Frozen scope or parameters changed')
    with mp.workdps(100):
        reviews=dict(exact_radius_and_scale_powers=exact_scale_review(),all_later_domain_threshold=threshold_review(),
                     global_radius_phase_and_frame=frame_review(),full_Haar_primary_covariance=covariance_review(),
                     independent_sourced_projected_ODE=sourced_review(),exact_higher_harmonic_factor=harmonic_factor_review(),
                     liftwise_mixed_band_sources=lifted_source_review(),fixed_order_jets_and_flat_tails=derivative_and_tail_review(),
                     failed_shortcuts=failed_shortcuts_review())
        refinement=refinement_review(reviews['full_Haar_primary_covariance'],covariance_review(True))
    inherited=inherited_provenance()
    gates={name:all(reviews[name]['checks'].values()) for name in REVIEW_NAMES}
    gates['refined_ODE_and_full_interval_quadrature_agree']=all(refinement['checks'].values())
    gates['all_inherited_bytes_preserved']=bool(inherited) and all(inherited.values())
    gates['conditional_scope_not_promoted']=not any(FLAGS[k] for k in
        ('actual_global_leading_profile_or_summed_background_numerically_evaluated',
         'actual_profile_value_gap_and_jet_constants_interval_certified',
         'leading_profile_background_and_mean_inverse_foundation_independently_verified',
         'full_global_nonlinear_hierarchy_independently_revalidated',
         'unconditional_complete_residual_improvement_verified','unconditional_blowup_construction_verified','blowup_verified'))
    return encode(dict(status=STATUS,protocol=expected_protocol(),flags=FLAGS,proof=proof_ledger(),
                       reviews=reviews,refinement=refinement,inherited_provenance=inherited,
                       source_hashes=source_hashes(),gates=gates,remaining_obligations=REMAINING))


def validate(record):
    if record.get('status')!=STATUS or record.get('flags')!=FLAGS or record.get('proof')!=proof_ledger():
        raise ValueError('Conditional pulse replay does not independently verify the global profile, residual or blow-up')
    if record.get('protocol')!=expected_protocol() or record.get('remaining_obligations')!=list(REMAINING):
        raise ValueError('Frozen scope or remaining obligations changed')
    if set(record.get('gates',{}))!=set(GATE_NAMES) or not all(record['gates'].values()):
        raise ValueError('A required estimate, diagnostic, scope or provenance gate failed')


def compare_records(actual,frozen,location='root'):
    if isinstance(actual,dict):
        if not isinstance(frozen,dict) or set(actual)!=set(frozen):
            raise ValueError('Frozen fields differ at '+location)
        for key in actual:
            compare_records(actual[key],frozen[key],location+'.'+key)
    elif isinstance(actual,list):
        if not isinstance(frozen,list) or len(actual)!=len(frozen):
            raise ValueError('Frozen sequence differs at '+location)
        for j,value in enumerate(actual):
            compare_records(value,frozen[j],location+'.'+str(j))
    elif location.startswith(('root.reviews.','root.refinement.')) and isinstance(actual,float):
        # Tiny full-Haar masses are about 1e-16: a blanket 1e-10 absolute
        # tolerance would erase their normalization. Roundoff-error counters
        # instead need an absolute tolerance across BLAS/SciPy versions.
        strict_scale=(location.endswith(('.h','.covariance_lower','.covariance_upper')) or '.column.' in location)
        absolute=1e-300 if strict_scale else 2e-10
        if isinstance(frozen,bool) or not isinstance(frozen,(float,int)) or not math.isclose(actual,frozen,rel_tol=3e-7,abs_tol=absolute):
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
    print(str(len(record['gates']))+' gates pass; '+str(len(record['inherited_provenance']))+' inherited bytes checked.')
    print('Actual leading/profile/background constants and the unconditional residual/blow-up construction remain unverified.')


if __name__=='__main__':
    main()
