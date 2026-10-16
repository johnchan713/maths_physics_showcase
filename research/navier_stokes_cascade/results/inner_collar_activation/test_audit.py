#!/usr/bin/env python3
"""Scoped activation, complete memory, endpoint and evidence regression tests."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from audit import HERE, compare_records, generate, main, validate
from collar import Activation, Jet, collar_threshold, exprel, flat_factor, flat_quotient, ga, sigma


class InnerActivationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.record = generate()

    def test_positive_quotient_matches_independent_finite_interval(self):
        self.assertTrue(all(self.record['reviews']['flat_quotient_and_endpoint_jets']['checks'].values()))

    def test_quotient_endpoint_is_coefficient_over_two_not_zero(self):
        for power in (0, 3, 6):
            self.assertEqual(flat_quotient(0, lambda u: 7+u, power=power), mp.mpf('3.5'))

    def test_exact_first_two_endpoint_jets(self):
        r = self.record['reviews']['flat_quotient_and_endpoint_jets']
        self.assertLess(r['endpoint_second_jet_scaled_error'], .001)

    def test_radial_jacobian_is_not_optional(self):
        self.assertGreater(self.record['reviews']['flat_quotient_and_endpoint_jets']['omitted_Jacobian_gap'], .01)

    def test_all_five_reference_primitive_and_eta_derivatives(self):
        self.assertTrue(all(self.record['reviews']['five_forward_primitive_derivatives']['checks'].values()))

    def test_complete_forward_map_difference_matches_literal_subtraction(self):
        r = self.record['reviews']['complete_activation_difference_and_eta_jets']
        self.assertTrue(r['checks']['complete_moment_map_matches_resolvable_subtraction'])

    def test_eta_jets_of_every_moment_difference(self):
        self.assertLess(self.record['reviews']['complete_activation_difference_and_eta_jets']['eta_jet_central_difference_error'], 1e-10)

    def test_pressure_and_angular_moment_memory_cannot_be_dropped(self):
        r = self.record['reviews']['complete_activation_difference_and_eta_jets']['rows'][1]
        self.assertGreater(abs(r['pressure_p2_hat']), 1e-14)
        self.assertGreater(abs(r['angular_Jeta_p1_hat']), 1e-14)

    def test_normalized_values_survive_raw_flat_underflow(self):
        self.assertEqual(float(flat_factor('.00001')), 0)
        self.assertGreater(flat_quotient('.00001', lambda u: 1+u), 0)
        r = self.record['reviews']['complete_activation_difference_and_eta_jets']['rows'][-1]
        self.assertGreater(abs(r['moment_coefficients'][0]), .01)

    def test_two_integrations_give_the_sixth_order_endpoint(self):
        self.assertTrue(all(self.record['reviews']['endpoint_coefficients_and_explicit_collar']['checks'].values()))

    def test_squared_cone_loss_uses_width_squared_and_sixth_power(self):
        rows = self.record['reviews']['complete_activation_difference_and_eta_jets']['rows']
        self.assertTrue(all(0 <= r['excess_cone_loss_over_t2s6'] < 10 for r in rows))

    def test_kappa_zero_endpoint_requires_no_inverse_kappa(self):
        r = Activation('.001', 0, order=16).stress('.05', 0)
        self.assertGreater(r['denominator_over_ea'], 0)
        self.assertTrue(mp.isfinite(r['cross_over_ea']))

    def test_entire_exprel_and_derivative_at_zero(self):
        z = exprel(Jet(0, 1))
        self.assertEqual(z.v, 1)
        self.assertEqual(z.d, mp.mpf('.5'))
        small = mp.mpf('1e-30')
        self.assertLess(abs(exprel(small).v-mp.expm1(small)/small), mp.mpf('1e-75'))

    def test_explicit_collar_budgets(self):
        r = collar_threshold('2.3', '10', '1000000')
        self.assertLessEqual(r['smax'], mp.mpf('.25'))
        self.assertLessEqual(r['error_bound']*r['width_max']*r['smax']**3, r['allowance']*(1+mp.mpf('1e-60')))
        self.assertLessEqual(sigma(r['smax'])*11, mp.mpf('.15'))

    def test_tiny_core_tolerance_is_not_exact_flat_relative_control(self):
        r = self.record['reviews']['baseline_defect_and_failed_shortcuts']
        self.assertGreater(r['log10_tiny_baseline_error_over_activation'], 900)
        self.assertGreater(r['erased_defect_gap'], 1)

    def test_actual_fixture_stress_keeps_reference_defect(self):
        self.assertTrue(self.record['reviews']['baseline_defect_and_failed_shortcuts']['checks']['complete_seed_defect_identity'])
        r = self.record['reviews']['complete_activation_difference_and_eta_jets']
        self.assertFalse(r['seed_is_actual_analytic_solution'])
        self.assertTrue(r['cone_test_is_activated_excess_only'])

    def test_growth_keeps_axial_shear_denominator(self):
        self.assertEqual(self.record['reviews']['baseline_defect_and_failed_shortcuts']['wrong_growth_factor'], 4)

    def test_refined_positive_integrals_agree(self):
        self.assertTrue(all(self.record['reviews']['tighter_positive_quadrature']['checks'].values()))

    def test_rescaled_ga_bounds(self):
        for s in (0, '.001', '.05', '.25'):
            self.assertGreaterEqual(ga(s), 2)
            self.assertLessEqual(ga(s), 6)

    def test_invalid_domains_and_false_seed_gaps_rejected(self):
        for width, kappa in ((0, 0), ('.02', 0), ('.01', '-.1'), ('.01', '.6')):
            with self.assertRaises(ValueError):
                Activation(width, kappa)
        for s in ('-.01', '.3'):
            with self.assertRaises(ValueError):
                flat_quotient(s, lambda u: 1)
        with self.assertRaises(ValueError):
            flat_quotient('.1', lambda u: 1, power=-1)
        with self.assertRaises(ValueError):
            collar_threshold(2, 10, 1)

    def test_inherited_hashes_and_frozen_evidence_reproduce(self):
        validate(self.record)
        self.assertTrue(all(self.record['inherited_provenance'].values()))
        compare_records(self.record, json.loads((HERE/'evidence.json').read_text()))

    def test_actual_seed_or_blowup_promotion_is_rejected(self):
        for key in ('actual_core_contraction_independently_reverified',
                    'complete_inner_and_compact_foundation_verified',
                    'unconditional_uniform_actual_pulse_estimate_verified',
                    'unconditional_complete_residual_improvement_verified', 'blowup_verified'):
            altered = deepcopy(self.record)
            altered['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(altered)

    def test_omitted_obligation_or_erased_coefficient_is_rejected(self):
        altered = deepcopy(self.record)
        altered['remaining_obligations'].pop()
        with self.assertRaises(ValueError):
            validate(altered)
        altered = deepcopy(self.record)
        altered['reviews']['complete_activation_difference_and_eta_jets']['rows'][0]['pressure_p2_hat'] = 0.
        with self.assertRaises(ValueError):
            compare_records(self.record, altered)

    def test_cli_never_overwrites_a_record(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'record.json'
            path.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output', str(path)])
            self.assertEqual(path.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
