#!/usr/bin/env python3
"""Current force transfer, exact localization, gluing, energy and scope controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from forcing import TERMS,requirements,first_controlled_stage,choose_log2_scale,extension_review
from current_checks import localization_stats,forced_pressure_stats
from audit import HERE,FLAGS,generate,validate,compare_records,main,historical_review,parent_review


class CurrentForceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def review(self,name):return self.record['reviews'][name]

    def historical(self,name):return historical_review()['controls'][name]['checks']

    def test_current_actual_endpoint_and_realization_are_input(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'e4e38b2309d3a9d426c820b659df2a1bf80e15fc')
        self.assertTrue(self.record['gates']['current_actual_endpoint_and_realization_replayed'])

    def test_historical_old_field_not_substituted_for_current_one(self):
        self.assertTrue(historical_review()['old_force_field_not_substituted_for_current'])
        self.assertTrue(all(historical_review()['checks'].values()))

    def test_strict_current_domain_margin_and_first_band_preserved(self):
        r=self.review('current_domain_and_support')
        self.assertEqual(Q(r['q_over_q_e']),Q(3,8))
        self.assertEqual(Q(r['fourQ_over_qstar']),Q(3,16))
        self.assertTrue(all(r['checks'].values()))

    def test_full_current_pressure_tails_and_differences_cut_without_new_gauge(self):
        self.assertTrue(self.review('current_domain_and_support')['checks']['all_base_pressure_tails_and_every_current_pressure_difference_cut_in_same_gauge'])

    def test_twelve_localization_terms_all_present(self):
        self.assertEqual(len(TERMS),12)
        self.assertEqual(tuple(self.review('complete_localization_and_pressure')['retained_terms']),TERMS)

    def test_complete_cut_velocity_product_rule_through_two_derivatives(self):
        self.assertLess(localization_stats()['velocity_product_error'],2e-10)

    def test_complete_current_q_localization_residual_identity(self):
        r=localization_stats();self.assertEqual(r['points'],4)
        self.assertEqual(r['positive_viscosities'],[.6,1.,1.7]);self.assertLess(r['identity_error'],2e-10)

    def test_complete_localized_velocity_divergence_is_zero(self):
        self.assertLess(localization_stats()['divergence_error'],2e-12)

    def test_naive_velocity_cutoff_breaks_divergence(self):
        self.assertGreater(localization_stats()['naive_divergence'],1e-4)

    def test_theta_dependent_cutoff_is_not_allowed_for_direct_swirl(self):
        self.assertGreater(localization_stats()['nonaxisymmetric_direct_swirl_divergence'],1e-4)

    def test_every_retained_term_has_nonzero_omission_control(self):
        r=localization_stats()['retained_term_sizes']
        self.assertEqual(tuple(r),TERMS);self.assertGreater(min(r.values()),1e-7)

    def test_complete_pressure_cutoff_gradient_cannot_be_discarded(self):
        self.assertGreater(localization_stats()['retained_term_sizes']['pressure_cutoff'],1e-7)

    def test_cutoff_time_and_all_viscous_rows_are_nonzero(self):
        r=localization_stats()['retained_term_sizes']
        for n in ('cutoff_time','curl_time','cutoff_laplacian','gradient_viscosity','curl_viscosity'):
            self.assertGreater(r[n],1e-7)

    def test_all_nonlinear_localization_rows_are_retained(self):
        r=localization_stats()['retained_term_sizes']
        for n in ('quadratic_rescaling','velocity_cutoff_transport','old_curl_transport','curl_old_transport','curl_square'):
            self.assertGreater(r[n],1e-7)

    def test_only_excluded_transport_vanishes_by_cross_product_identity(self):
        self.assertLess(localization_stats()['orthogonality_error'],2e-12)

    def test_independent_polynomial_forced_Poisson_identity(self):
        self.assertLess(forced_pressure_stats()['forced_Poisson_error'],2e-12)

    def test_zero_velocity_divergence_does_not_mean_zero_force_divergence(self):
        r=forced_pressure_stats();self.assertLess(abs(r['velocity_divergence']),2e-12)
        self.assertGreater(abs(r['force_divergence']),1e-3);self.assertGreater(r['omitted_divergence_gap'],1e-3)

    def test_all_fixed_force_jet_requirements_and_next_time_order(self):
        for m in (0,5,1000):
            r=requirements(m)
            self.assertEqual(r['potential_bound'],m+3);self.assertEqual(r['direct_swirl_bound'],m+2)
            self.assertEqual(r['pressure_bound'],m+1);self.assertEqual(r['next_time_potential_bound'],m+4)

    def test_invalid_derivative_orders_rejected(self):
        for m in (-1,True,1.5):
            with self.assertRaises(ValueError):requirements(m)

    def test_complete_residual_corner_flatness_not_only_regular_angular_source(self):
        c=self.review('global_terminal_force_jets')['checks']
        self.assertTrue(c['complete_residual_not_only_regular_angular_part_is_flat'])
        self.assertTrue(c['compact_X_uniformity_plus_exact_heat_exterior_zero_used'])
        self.assertTrue(c['origin_flatness_alone_not_used_for_away_from_origin_limits'])

    def test_first_controlled_stage_includes_scalar_case(self):
        self.assertEqual(first_controlled_stage(0,0),1)
        self.assertEqual(first_controlled_stage(1,3),8)

    def test_extreme_synthetic_extension_constraints_met(self):
        r=extension_review();self.assertTrue(all(r['checks'].values()))
        self.assertGreater(r['rows'][-1]['dyadic_exponent'],100000)
        self.assertTrue(r['actual_terminal_norms_and_scales_not_evaluated'])

    def test_invalid_extension_stage_and_previous_power_rejected(self):
        for args in ((0,[(0,0,0)]),(True,[(0,0,0)]),(1,[],1),(1,[(0,0,0)],0)):
            with self.assertRaises(ValueError):choose_log2_scale(*args)

    def test_invalid_extension_mixed_orders_and_constants_rejected(self):
        for cs in ([(0,1,0)],[(0,0,-1)],[(0,0,float('inf'))],[(True,0,0)]):
            with self.assertRaises(ValueError):choose_log2_scale(1,cs)

    def test_every_fixed_mixed_tail_eventually_summable(self):
        for a,m in ((0,0),(20,4),(6,70)):
            J=first_controlled_stage(a,m)
            for j in (J,J+1,J+100):
                self.assertLessEqual(a+m,j//2);self.assertLess(m,j)

    def test_historical_mixed_terminal_jet_gluing_and_compact_support_reproduced(self):
        c=self.historical('terminal_mixed_jet_gluing')
        self.assertTrue(c['all_available_mixed_terminal_jets_match'])
        self.assertTrue(c['spatial_support_time_support_and_flat_origin_are_preserved'])

    def test_zero_future_force_shortcut_fails(self):
        self.assertTrue(self.historical('terminal_mixed_jet_gluing')['extending_the_force_by_zero_would_create_a_jump'])

    def test_failed_low_energy_quadrature_is_preserved_not_tolerance_relaxed(self):
        c=self.historical('whole_space_energy_and_dissipation')
        self.assertTrue(c['low_quadrature_failure_is_retained'])
        self.assertTrue(c['independent_viscous_integration_by_parts'])

    def test_independent_complete_forced_energy_and_flux_controls_reproduced(self):
        c=self.historical('whole_space_energy_and_dissipation')
        for n in ('complete_forced_energy_identity','pressure_and_nonlinear_integrals_cancel',
            'integrated_energy_bound_holds','omission_of_dissipation_is_detected'):self.assertTrue(c[n])

    def test_zero_norm_regularization_and_monotone_dissipation_used(self):
        c=self.review('current_whole_space_energy_and_core_growth')['checks']
        self.assertTrue(c['no_division_by_possibly_zero_velocity_norm'])
        self.assertTrue(c['monotone_convergence_only_applied_to_nonnegative_dissipation'])

    def test_energy_bound_retains_factor_two_and_only_preterminal_solution(self):
        r=self.review('current_whole_space_energy_and_core_growth')
        self.assertIn('+2*integral_0^t',r['bound'])
        self.assertTrue(r['checks']['no_global_competitor_comparison_or_regular_velocity_at_t_one_inferred'])

    def test_official_admissibility_requires_full_decay_not_zero_force_divergence(self):
        r=self.review('official_admissibility_interface')
        self.assertTrue(r['checks']['no_divergence_free_force_requirement_added_to_official_definition'])
        self.assertTrue(r['checks']['compact_smooth_force_implies_every_required_spacetime_decay_bound'])

    def test_smooth_future_force_does_not_construct_future_velocity_or_Taylor_convergence(self):
        self.assertTrue(FLAGS['conditional_one_smooth_compact_global_force_matching_exact_preterminal_residual_supplied'])
        self.assertFalse(FLAGS['smooth_velocity_or_pressure_after_t_one_supplied'])
        self.assertFalse(FLAGS['raw_force_Taylor_convergence_or_uniform_terminal_order_constants_assumed'])

    def test_current_core_and_Gaussian_exclusion_preserved(self):
        self.assertTrue(parent_review()[0]['flags']['actual_flat_velocities_base_pressure_tails_and_Gaussian_exclusion_preserved'])
        self.assertTrue(FLAGS['conditional_complete_core_lower_bound_survives_current_localization'])
        self.assertFalse(FLAGS['Gaussian_cutoff_tails_forward_inverted'])

    def test_732_inherited_hashes_and_16_own_source_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),732)
        self.assertTrue(all(self.record['inherited_provenance'].values()))
        self.assertEqual(len(self.record['source_hashes']),16)

    def test_frozen_record_reproduces(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r)

    def test_record_roundtrips_without_numpy_booleans_or_roundoff_refreezing(self):
        r=json.loads(json.dumps(self.record));validate(r);compare_records(self.record,r)
        self.assertTrue(self.review('current_q_localization_and_forced_pressure_controls')['raw_roundoff_scores_tested_not_refrozen'])

    def test_unconditional_competitor_exclusion_and_blowup_promotions_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError):validate(r)

    def test_analytic_foundation_premises_cannot_be_erased(self):
        r=deepcopy(self.record)
        r['flags']['inherited_leading_background_inverse_and_uniform_physical_estimates_remain_analytical_premises']=False
        with self.assertRaises(ValueError):validate(r)

    def test_mutated_pressure_or_energy_ledger_rejected_despite_passed_gates(self):
        for n in ('pressure','energy'):
            r=deepcopy(self.record)
            if n=='pressure':r['reviews']['complete_localization_and_pressure']['pressure_identity']='-Delta p=sum dd(uu)'
            else:r['reviews']['current_whole_space_energy_and_core_growth']['bound']='||u||^2+integral grad^2<=F^2'
            with self.assertRaises(ValueError):validate(r)

    def test_missing_gate_obligation_or_provenance_rejected(self):
        for section in ('gates','remaining_obligations','source_hashes','inherited_provenance'):
            r=deepcopy(self.record)
            if section=='remaining_obligations':r[section].pop()
            else:r[section].pop(next(iter(r[section])))
            with self.assertRaises(ValueError):validate(r)

    def test_cli_does_not_overwrite_frozen_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':unittest.main()
