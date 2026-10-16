#!/usr/bin/env python3
"""Phase-locked covariance, fixed signed inverse, pressure and scope checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE, FLAGS, generate, validate, compare_records, main
from signed import KAPPA, fixture, orders, physical_powers


class RealSignedUpdateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = generate()

    def review(self, name):
        return self.record['reviews'][name]

    def test_full_tensor_increment_and_independent_angular_moments(self):
        r = self.review('exact_phase_locked_full_curl_covariance')
        self.assertEqual(r['exact_tensor_entry_checks'],36)
        self.assertTrue(r['checks']['exact_full_tensor_increment'])
        self.assertTrue(r['checks']['independent_quarter_phase_tensor_average'])

    def test_old_and_new_complete_curls_are_divergence_free(self):
        r = self.review('exact_phase_locked_full_curl_covariance')
        self.assertEqual(r['exact_divergence_checks'],12)
        self.assertTrue(r['checks']['primary_and_signed_complete_curl_divergence'])

    def test_all_six_updated_covariances_are_retained(self):
        self.assertTrue(self.review('exact_phase_locked_full_curl_covariance')['checks']['all_six_updated_covariances_nonzero'])

    def test_variable_amplitude_curl_is_not_amplitude_times_old_curl(self):
        self.assertTrue(self.review('exact_phase_locked_full_curl_covariance')['checks']['signed_amplitude_curl_commutators_retained'])

    def test_both_new_square_and_old_new_curl_cross_are_nonzero(self):
        checks = self.review('exact_phase_locked_full_curl_covariance')['checks']
        self.assertTrue(checks['signed_increment_square_not_discarded'])
        self.assertTrue(checks['primary_new_sine_curl_cross_not_discarded'])
        self.assertTrue(checks['total_pair_remainder_is_not_zero'])

    def test_only_two_principal_cross_entries_are_prescribed(self):
        _,_,delta,principal,remainder = fixture()
        self.assertNotEqual(principal[2][2],0)
        self.assertEqual(delta[2][2],principal[2][2]+remainder[2][2])

    def test_general_old_cosine_quadrature_destroys_special_cancellation(self):
        self.assertTrue(self.review('exact_phase_locked_full_curl_covariance')['checks']['generic_old_cosine_error_breaks_sharp_identity'])
        self.assertFalse(FLAGS['sharper_covariance_formula_applies_to_arbitrary_later_complex_amplitudes'])

    def test_fixture_orthogonality_does_not_replace_actual_label_separation(self):
        self.assertTrue(self.review('exact_phase_locked_full_curl_covariance')['distinct_fixture_frequencies_not_replacement_for_actual_closed_support_separation'])

    def test_signed_inverse_has_negative_coefficient_and_no_sign_requirement(self):
        r = self.review('signed_inverse_edge_and_partition_controls')
        self.assertEqual(r['inverse_coefficients'],['-1450/7','100/7'])
        self.assertTrue(r['checks']['fixed_signed_linear_inverse_recovers_target'])

    def test_real_cosine_factor_two_is_required(self):
        self.assertTrue(self.review('signed_inverse_edge_and_partition_controls')['checks']['factor_two_is_necessary'])

    def test_physical_partition_is_squared_not_divided_out(self):
        r = self.review('signed_inverse_edge_and_partition_controls')['checks']
        self.assertTrue(r['squared_partition_physical_normalization'])
        self.assertTrue(r['no_division_by_slow_partition_needed'])

    def test_reciprocal_root_derivative_half_weight_calculus(self):
        r = self.review('signed_inverse_edge_and_partition_controls')
        self.assertEqual(r['exact_edge_term_checks'],13)
        self.assertTrue(r['checks']['every_checked_reciprocal_root_term_retains_half_edge_weight'])
        for n in range(100):
            self.assertEqual((-n-Q(1,2))+n+1,Q(1,2))

    def test_weighted_bounds_do_not_assert_positive_updated_amplitudes(self):
        self.assertTrue(self.review('signed_inverse_edge_and_partition_controls')['checks']['weighted_class_does_not_ensure_positive_updated_amplitudes'])

    def test_recomputed_pressure_compact_radial_identity_and_moment(self):
        r = self.review('recomputed_axial_pressure_and_fluxes')['checks']
        self.assertTrue(all(r[n] for n in ('delta_pressure_compact_at_both_ends','delta_pressure_radial_source','whole_delta_axial_pressure_moment')))

    def test_both_integrated_updated_flux_identities(self):
        r = self.review('recomputed_axial_pressure_and_fluxes')['checks']
        self.assertTrue(r['delta_angular_integrated_balance'])
        self.assertTrue(r['delta_axial_integrated_balance'])

    def test_stress_cancellation_retains_axial_flux_and_pressure(self):
        r = self.review('recomputed_axial_pressure_and_fluxes')['checks']
        self.assertTrue(r['angular_stress_cancellation_exact'])
        self.assertTrue(r['axial_stress_cancellation_exact'])
        self.assertTrue(r['omitted_delta_axial_pressure_detected'])
        self.assertTrue(r['omitted_delta_axial_covariance_detected'])

    def test_pressure_bump_moves_and_all_defects_change(self):
        r = self.review('recomputed_axial_pressure_and_fluxes')['checks']
        self.assertTrue(r['moving_bump_product_derivative_nonzero'])
        self.assertTrue(r['covariance_update_changes_all_three_defects'])

    def test_initial_covariance_gain_exact_and_limited_by_signed_square(self):
        r = orders(2-3*KAPPA)
        self.assertEqual(r['real_pair_remainder'],3-6*KAPPA)
        self.assertEqual(r['signed_transverse_square'],r['real_pair_remainder'])
        self.assertEqual(r['real_covariance_gain'],1-3*KAPPA)

    def test_alpha_one_square_blocks_a_positive_gain(self):
        self.assertEqual(orders(1)['real_covariance_gain'],0)
        with self.assertRaises(ValueError):
            orders(Q(9,10))

    def test_covariance_gain_formula_and_minimum_for_both_regimes(self):
        for alpha in (Q(11,10),Q(3,2),Q(2),Q(3),Q(20)):
            self.assertEqual(orders(alpha)['real_covariance_gain'],min(1-2*KAPPA,alpha-1))

    def test_sharp_average_does_not_improve_full_mean_or_preserve_S2(self):
        r = orders(2-3*KAPPA)
        self.assertEqual(r['full_tangential_mean'],2-4*KAPPA)
        self.assertEqual(r['recomputed_defects'],2-3*KAPPA)
        self.assertFalse(FLAGS['complete_unaveraged_mean_residual_already_improved'])
        self.assertFalse(FLAGS['compatibility_defects_asserted_to_remain_in_S2_after_wave_change'])

    def test_changed_pressure_has_distinct_full_and_averaged_orders(self):
        r = orders(2-3*KAPPA)
        self.assertEqual(r['averaged_pressure_change'],2-3*KAPPA)
        self.assertEqual(r['full_pressure_change'],2-4*KAPPA)
        self.assertEqual(r['axial_flux_and_averaged_pressure'],3-3*KAPPA)

    def test_this_substep_generated_sources_include_viscosity_and_all_interactions(self):
        r = orders(2-3*KAPPA)
        self.assertEqual(r['linear_wave_remainder'],2-6*KAPPA)
        self.assertEqual(r['actual_mean_cross'],2-4*KAPPA)
        self.assertEqual(r['signed_wave_self_interaction'],3-7*KAPPA)

    def test_uniform_signed_physical_normalizations(self):
        alpha = 2-3*KAPPA; h = Q(1,200)
        p = physical_powers(alpha,h,0,0)
        self.assertEqual(p,dict(velocity=Q(-9950003,20000000),pressure=Q(-20000003,20000000),potential=Q(99997,20000000)))
        for a,b in ((1,0),(0,1),(2,3),(6,4)):
            shift = -a*Q(201,400)-b*Q(403,400)
            self.assertEqual(physical_powers(alpha,h,a,b),{key:value+shift for key,value in p.items()})

    def test_invalid_fixed_parameters_and_derivative_orders_rejected(self):
        for k in (0,-1,Q(1,2)):
            with self.assertRaises(ValueError):
                orders(2,k)
        for h,a,b in ((0,0,0),(Q(1,100),0,0),(Q(1,200),-1,0),(Q(1,200),True,0)):
            with self.assertRaises(ValueError):
                physical_powers(2,h,a,b)

    def test_normalized_covariance_gain_is_not_absolute_physical_decay(self):
        nu = 3-6*KAPPA; h = Q(1,200); A = Q(1,2)+h
        self.assertEqual(-2*A+h*nu,Q(-9950003,10000000))
        self.assertLess(-2*A+h*nu,0)

    def test_actual_parent_and_610_historical_hashes_replayed(self):
        self.assertTrue(self.record['gates']['actual_initial_means_parent_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),610)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_record_reproduces(self):
        frozen = json.loads((HERE/'evidence.json').read_text())
        validate(frozen);compare_records(self.record,frozen)

    def test_general_cycle_full_residual_and_blowup_promotions_rejected(self):
        for name,value in FLAGS.items():
            if not value:
                r = deepcopy(self.record);r['flags'][name] = True
                with self.assertRaises(ValueError):
                    validate(r)

    def test_phase_locking_and_analytical_dependencies_cannot_be_removed(self):
        for section,name in (('protocol','same_phase_in_each_old_new_label'),
                             ('flags','relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas')):
            r = deepcopy(self.record);r[section][name] = False
            with self.assertRaises(ValueError):
                validate(r)

    def test_missing_gate_or_obligation_rejected(self):
        for section in ('gates','remaining_obligations'):
            r = deepcopy(self.record)
            if section == 'gates':
                r[section].pop(next(iter(r[section])))
            else:
                r[section].pop()
            with self.assertRaises(ValueError):
                validate(r)

    def test_changed_diagnostic_or_exponent_rejected(self):
        for review,section,key,value in (
            ('signed_inverse_edge_and_partition_controls','inverse_coefficients',0,'0'),
            ('initial_covariance_source_and_physical_orders','rows','real_pair_remainder','3')):
            r = deepcopy(self.record);r['reviews'][review][section][key] = value
            with self.assertRaises(ValueError):
                validate(r)

    def test_changed_source_or_provenance_rejected(self):
        for section in ('source_hashes','inherited_provenance'):
            r = deepcopy(self.record);r[section][next(iter(r[section]))] = False
            with self.assertRaises(ValueError):
                validate(r)

    def test_cli_never_overwrites_frozen_record(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__ == '__main__':
    unittest.main()
