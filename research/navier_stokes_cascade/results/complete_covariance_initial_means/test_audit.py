#!/usr/bin/env python3
"""Complete covariance, pressure memory, initial Haar gain and scope controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE, FLAGS, compare_records, generate, main, validate
from means import Poly, R, Z, U, fixture, raw_data, KAPPA


class InitialMeanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = generate()

    def review(self, name):
        return self.record['reviews'][name]

    def test_exact_arithmetic_rejects_float_inputs(self):
        with self.assertRaises(TypeError):
            Poly(.1)
        self.assertEqual((R+Q(1,3))**2, R**2+Q(2,3)*R+Q(1,9))

    def test_polynomial_primitive_and_endpoints(self):
        p = (R-1)**3*(2-R)**3*(1+Z)
        self.assertEqual(p.primitive().diff(0), p)
        self.assertEqual(p.at(0,1), 0)
        self.assertEqual(p.at(0,2), 0)

    def test_unsupported_logarithmic_primitive_rejected(self):
        with self.assertRaises(ValueError):
            Poly(1).shift(-1).integral()

    def test_auxiliary_haar_moments_are_exact_not_quadrature(self):
        self.assertEqual(U.haar(), 0)
        self.assertEqual((U**2).haar(), Q(1,2))
        self.assertEqual(((R+Z*U)**2).haar(), R**2+Z**2/2)
        with self.assertRaises(ValueError):
            (U**3).haar()

    def test_complete_curl_divergence_and_six_covariances(self):
        checks = self.review('exact_complete_curl_and_pressure_moments')['checks']
        self.assertTrue(all(checks[n] for n in ('complete_cosine_divergence_zero',
            'complete_sine_divergence_zero','covariance_symmetric','all_six_covariances_nonzero')))

    def test_actual_fixture_does_not_drop_axial_velocity(self):
        _, _, c, w = fixture()
        self.assertNotEqual(c[2], 0)
        self.assertNotEqual(w[2][2], 0)
        self.assertNotEqual(w[1][2], 0)
        self.assertTrue(self.review('exact_complete_curl_and_pressure_moments')['checks']['omitted_axial_velocity_detected'])

    def test_raw_defect_integrations_by_parts(self):
        checks = self.review('exact_complete_curl_and_pressure_moments')['checks']
        self.assertTrue(checks['raw_pressure_defect_ibp'])
        self.assertTrue(checks['raw_axial_defect_ibp'])

    def test_exact_nonzero_rational_pressure_control(self):
        r = self.review('exact_complete_curl_and_pressure_moments')
        self.assertEqual(r['nonzero_controls']['P_at_Z0'], '201/2974571600')
        self.assertEqual(r['nonzero_controls']['moving_bump_cprime'], '1/40')

    def test_compact_pressure_zero_at_both_endpoints_and_radial_identity(self):
        r = self.review('exact_complete_curl_and_pressure_moments')['checks']
        self.assertTrue(r['moving_bump_mass_one'])
        self.assertTrue(r['compact_pressure_at_both_ends'])
        self.assertTrue(r['compact_pressure_radial_identity'])

    def test_axial_pressure_moment_and_moving_product_derivative(self):
        r = self.review('exact_complete_curl_and_pressure_moments')['checks']
        self.assertTrue(r['compact_pressure_axial_moment'])
        self.assertTrue(r['whole_moving_bump_derivative'])
        self.assertTrue(r['omitted_bump_derivative_detected'])
        self.assertTrue(r['omitted_axial_pressure_detected'])

    def test_haar_zero_means_remove_fast_slow_cross_terms_only_after_averaging(self):
        r = self.review('exact_Haar_and_recomputed_defect_changes')['checks']
        self.assertTrue(r['all_fast_haar_means_zero'])
        self.assertTrue(r['old_fast_slow_products_have_zero_mean'])
        self.assertNotEqual(R*U*Z, 0)
        self.assertEqual((R*U*Z).haar(), 0)

    def test_fast_linear_changes_cancel_but_fast_squares_remain(self):
        self.assertTrue(self.review('exact_Haar_and_recomputed_defect_changes')['checks']['fast_averaged_linear_changes_zero'])
        self.assertNotEqual(((R*U)**2).haar(), 0)

    def test_compact_slow_curl_and_both_velocity_moments(self):
        r = self.review('exact_Haar_and_recomputed_defect_changes')['checks']
        self.assertTrue(r['exact_slow_divergence'])
        self.assertTrue(r['slow_velocity_moments_zero'])

    def test_all_recomputed_defect_differences(self):
        r = self.review('exact_Haar_and_recomputed_defect_changes')['checks']
        self.assertTrue(all(r[n] for n in ('radial_recomputed_remainder','pressure_defect_difference',
                                          'angular_defect_difference','axial_defect_difference')))

    def test_slow_linear_map_does_not_set_nonlinear_defects_zero(self):
        self.assertTrue(self.review('exact_Haar_and_recomputed_defect_changes')['checks']['nonlinear_defects_not_identically_zero'])

    def test_sharp_initial_gain_not_valid_for_arbitrary_later_means(self):
        self.assertTrue(self.review('exact_Haar_and_recomputed_defect_changes')['checks']['nonzero_old_mean_breaks_sharp_reduction'])
        self.assertTrue(self.record['protocol']['initial_zero_Haar_cancellation_not_a_later_cycle_estimate'])

    def test_initial_raw_moments_have_no_graph_loss(self):
        rows = self.review('all_fixed_order_initial_exponent_ledger')['rows']
        self.assertEqual(Q(rows['raw_defects']), 1)
        self.assertEqual(Q(rows['raw_averaged_gr']), 1)
        self.assertEqual(Q(rows['raw_full_gr']), 1-KAPPA)

    def test_one_power_initial_compatibility_gain(self):
        rows = self.review('all_fixed_order_initial_exponent_ledger')['rows']
        self.assertEqual(Q(rows['final_three_defects'])-Q(rows['raw_defects']), 1)
        self.assertTrue(self.review('all_fixed_order_initial_exponent_ledger')['checks']['all_averaged_slow_terms_at_least_two'])

    def test_viscosity_graph_loss_not_erased_by_haar_moment_gain(self):
        rows = self.review('all_fixed_order_initial_exponent_ledger')['rows']
        self.assertEqual(Q(rows['fast_viscosity']), 2-3*KAPPA)
        self.assertEqual(Q(rows['final_full_tangential_residual']), 2-3*KAPPA)
        self.assertLess(Q(rows['final_full_tangential_residual']), 2)

    def test_pressure_tail_exception_and_cutoff_transition_retained(self):
        r = self.review('complete_background_pressure_tail_controls')
        self.assertTrue(all(r['checks'].values()))
        self.assertEqual(r['exact_controls']['second_tail_axial_numerator'], '-99/125')
        self.assertNotEqual(Q(r['exact_controls']['first_tail_transition']), 0)

    def test_old_maps_replayed_without_rewriting_historical_scope(self):
        self.assertTrue(self.record['gates']['retained_compact_mean_maps_replayed'])
        self.assertTrue(self.record['compact_map_replay']['inherited_conditional_scope_not_silently_changed'])

    def test_actual_primary_parent_and_595_historical_hashes(self):
        self.assertTrue(self.record['gates']['actual_primary_parent_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']), 595)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_evidence_reproduces(self):
        frozen = json.loads((HERE/'evidence.json').read_text())
        validate(frozen); compare_records(self.record, frozen)

    def test_later_covariance_source_residual_and_blowup_promotions_rejected(self):
        for name, value in FLAGS.items():
            if not value:
                r = deepcopy(self.record); r['flags'][name] = True
                with self.assertRaises(ValueError):
                    validate(r)

    def test_analytical_dependency_and_flat_scope_cannot_be_removed(self):
        for section, name in (('flags','relies_on_replayed_background_geometry_pulse_and_compact_mean_lemmas'),
                              ('protocol','retained_flat_base_residual_separate_from_weighted_shell_class')):
            r = deepcopy(self.record); r[section][name] = False
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

    def test_changed_evidence_controls_or_exponents_rejected(self):
        for review, section, key, value in (
            ('exact_complete_curl_and_pressure_moments','nonzero_controls','P_at_Z0','0'),
            ('all_fixed_order_initial_exponent_ledger','rows','final_full_tangential_residual','2')):
            r = deepcopy(self.record); r['reviews'][review][section][key] = value
            with self.assertRaises(ValueError):
                validate(r)

    def test_changed_current_or_historical_sources_rejected(self):
        for section in ('source_hashes','inherited_provenance'):
            r = deepcopy(self.record); r[section][next(iter(r[section]))] = False
            with self.assertRaises(ValueError):
                validate(r)

    def test_cli_never_overwrites_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'record.json'; path.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output', str(path)])
            self.assertEqual(path.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
