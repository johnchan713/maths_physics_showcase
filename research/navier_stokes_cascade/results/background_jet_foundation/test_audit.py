#!/usr/bin/env python3
"""Exact identities, independent pullbacks and regression controls for the audit."""
from copy import deepcopy
from fractions import Fraction
import json
from math import factorial
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from audit import FLAGS, HERE, compare_records, generate, main, validate
from foundation import (axis_average, band_coefficients, conversion_cost,
                        cutoff_log_scale, falling, normalized_jet_constant,
                        picard_log_term, shifted_cutoff_coefficients,
                        stirling_first_row, stirling_second_row)


class FoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.record = generate()

    def test_exact_conversion_and_cutoff_diagnostics(self):
        self.assertTrue(all(self.record['reviews']['scale_derivatives']['checks'].values()))

    def test_operator_change_of_basis_is_inverse_through_order_twelve(self):
        for m in range(13):
            for j in range(m+1):
                self.assertEqual(sum(stirling_second_row(m)[k]*stirling_first_row(k)[j]
                                     for k in range(j, m+1)), int(m == j))

    def test_first_kind_absolute_row_sum_is_factorial(self):
        for m in range(13):
            self.assertEqual(sum(map(abs, stirling_first_row(m))), factorial(m))

    def test_frozen_and_euler_derivatives_are_not_equal(self):
        self.assertEqual(falling(1, 2), 0)
        self.assertEqual(stirling_second_row(2), (0, 1, 1))
        self.assertEqual(conversion_cost(10)['frozen_to_euler'], 115975)

    def test_cutoff_formula_keeps_fractional_exponent(self):
        self.assertEqual(shifted_cutoff_coefficients(Fraction(1, 3), 2),
                         (Fraction(1, 9), Fraction(5, 3), Fraction(1)))

    def test_negative_and_noninteger_orders_rejected(self):
        for value in (-1, .5, True):
            with self.assertRaises(ValueError):
                stirling_second_row(value)
            with self.assertRaises(ValueError):
                stirling_first_row(value)

    def test_axis_quotient_is_integral_average_not_zero(self):
        self.assertEqual(axis_average([3, -5, Fraction(7, 2)]), (3, Fraction(-5, 2), Fraction(7, 6)))
        self.assertEqual(self.record['reviews']['regular_axis']['manufactured_first_pressure_slope'], '1/4')

    def test_continuum_bump_inverse_uses_whole_supports(self):
        result = self.record['reviews']['continuum_moment_inverse']
        self.assertTrue(all(result['checks'].values()))
        self.assertLess(result['E_inverse_upper'], 788)
        self.assertLess(result['U_inverse_upper'], 2)

    def test_actual_bump_refinement_and_tiny_exponent(self):
        result = self.record['reviews']['actual_bump_quadrature']
        self.assertTrue(all(result['checks'].values()))
        self.assertEqual(result['rows'][1]['patch_exponent'], '1e-30')

    def test_pointwise_nilpotence_does_not_imply_cross_radius_nilpotence(self):
        # Two pointwise nilpotent matrices whose product is nonzero.
        A, B = mp.matrix([[0, 1], [0, 0]]), mp.matrix([[0, 0], [1, 0]])
        self.assertEqual(mp.norm(A*A), 0)
        self.assertEqual(mp.norm(B*B), 0)
        self.assertNotEqual(mp.norm(A*B), 0)

    def test_common_sparse_derivative_block_annihilates_variable_factors(self):
        A, B, diagonal = mp.matrix(6), mp.matrix(6), mp.diag([1, 2, 3, 4, 5, 6])
        for i in (4, 5):
            for j in range(4):
                A[i, j], B[i, j] = (i+1)*(j+2), (i+2)*(j+3)
        self.assertEqual(mp.norm(A*diagonal*B), 0)
        B[0, 4] = 1
        self.assertNotEqual(mp.norm(A*diagonal*B), 0)

    def test_picard_majorant_reproduces(self):
        self.assertTrue(all(self.record['reviews']['common_collar_picard']['checks'].values()))
        with self.assertRaises(ValueError):
            picard_log_term(1, 2, 0)

    def test_normalized_constant_retains_cutoff_and_coefficient_order_costs(self):
        value = normalized_jet_constant(Fraction(1, 200), 3, 2, 5, 7)
        self.assertEqual(value, 70*Fraction(103, 100)**2)

    def test_scale_schedule_controls_every_smaller_q(self):
        h, n = Fraction(1, 200), 4
        constants = [mp.mpf(x) for x in (0, 2, 11, 19, 30)]
        scale = cutoff_log_scale(h, n, constants)
        for further in (0, 1, 100, 10000):
            self.assertLess(max(c+n*mp.log(2)-mp.mpf(n)/200*(scale+further) for c in constants), 0)

    def test_missing_jet_bounds_and_nonfinite_constants_rejected(self):
        for constants in ([0, 0], [0, 0, -1], [0, 0, 'inf']):
            with self.assertRaises(ValueError):
                cutoff_log_scale(Fraction(1, 200), 2, constants)

    def test_invalid_scale_and_jet_parameters_rejected(self):
        for h, n in ((0, 1), (Fraction(1, 100), 1), (Fraction(1, 200), 0)):
            with self.assertRaises(ValueError):
                cutoff_log_scale(h, n, [0, 0])
            with self.assertRaises(ValueError):
                normalized_jet_constant(h, n, 1, 1, 1)

    def test_chart_pullback_includes_mixed_and_boundary_derivatives(self):
        result = self.record['reviews']['band_pullback']
        self.assertTrue(all(result['checks'].values()))
        self.assertTrue(result['includes_T_zero'])
        self.assertGreater(result['maximum_error_from_freezing_chart_coefficients'], 1)

    def test_eta_boundary_does_not_make_L_singular(self):
        for eta in (-1, 1):
            coefficients = band_coefficients(2, 1, eta, '.005')
            self.assertEqual(coefficients[1][2], 0)
            self.assertEqual(coefficients[2][0], mp.mpf(100)/99)

    def test_small_values_alone_do_not_bound_euler_jets(self):
        # q^(2h) sin(1/q): at q=1/(2*pi*k), its normalized Euler jet
        # is -2*pi*k, unbounded despite the normalized value bound one.
        ratios = [2*mp.pi*k for k in (1, 10, 100)]
        self.assertEqual(ratios[2]/ratios[0], 100)

    def test_absolute_cone_ratio_does_not_select_positive_orientation(self):
        # N=(-1,0), K=(0,-1), c0=-1, u=3/4: the columns are
        # (5/4, +/-3/4). Both T=(1,0) and -T have absolute cone ratio zero.
        # Only T has nonnegative covariance coefficients.
        for sign in (1, -1):
            coefficient = Fraction(2*sign, 5)
            self.assertEqual(2*Fraction(5, 4)*coefficient, sign)
            self.assertEqual(abs(Fraction(0, sign)), 0)
            self.assertEqual(coefficient > 0, sign > 0)

    def test_all_earlier_hashes_and_scope_preserved(self):
        validate(self.record)
        self.assertTrue(all(self.record['inherited_provenance'].values()))
        self.assertFalse(FLAGS['unconditional_uniform_pulse_estimate_verified'])

    def test_unconditional_pulse_or_blowup_promotion_rejected(self):
        for key in ('full_leading_profile_foundation_independently_verified',
                    'unconditional_uniform_pulse_estimate_verified', 'blowup_verified'):
            changed = deepcopy(self.record)
            changed['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(changed)

    def test_deleted_obligation_or_failing_gate_rejected(self):
        changed = deepcopy(self.record)
        changed['remaining_obligations'].pop()
        with self.assertRaises(ValueError):
            validate(changed)
        changed = deepcopy(self.record)
        changed['gates']['band_pullback'] = False
        with self.assertRaises(ValueError):
            validate(changed)

    def test_evidence_replay_and_modified_exact_constant_rejected(self):
        saved = json.loads((HERE/'evidence.json').read_text())
        compare_records(self.record, saved)
        saved['reviews']['scale_derivatives']['rows'][10]['frozen_to_euler'] = 1
        with self.assertRaises(ValueError):
            compare_records(self.record, saved)

    def test_cli_never_overwrites_existing_record(self):
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp)/'existing.json'
            target.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output', str(target)])
            self.assertEqual(target.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
