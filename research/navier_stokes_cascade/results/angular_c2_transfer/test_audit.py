#!/usr/bin/env python3
"""Tests of C2 calculus, actual parameter inequalities, and limited promotion."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest
import mpmath as mp
from bounds import (HERE, core_transfer, normalization_constants, matching_transfer,
                    lower, upper, previous)
from jets import Jet2, physical_to_rows
from review import (normalization_review, cancellation_review, exponential_review,
                    axis_integrability_review)
from audit import STATUS, FLAGS, RESOLVED, REMAINING, validate, validate_protocol, inherited_provenance


class JetTests(unittest.TestCase):
    def test_product_retains_factor_two(self):
        e = Jet2(F(2, 3), F(1))
        self.assertEqual((e*e).dd, 2)

    def test_reciprocal_keeps_curvature(self):
        v = Jet2(F(3), F(2), F(5))
        self.assertEqual(v.inverse(), Jet2(F(1, 3), F(-2, 9), F(-7, 27)))

    def test_inverse_identity_is_exact(self):
        v = Jet2(F(3), F(2), F(5))
        self.assertEqual(v*v.inverse(), Jet2(F(1)))

    def test_zero_denominator_rejected(self):
        with self.assertRaises(ZeroDivisionError):
            Jet2(F(0)).inverse()

    def test_log_of_nonpositive_rejected(self):
        for v in (0, -1):
            with self.assertRaises(ValueError):
                Jet2(v).log()

    def test_unknown_input_derivatives_rejected(self):
        with self.assertRaises(ValueError):
            physical_to_rows([0]*5, 0, 16, 10000)

    def test_all_five_moments_required(self):
        with self.assertRaises(ValueError):
            physical_to_rows([Jet2(F(0))]*4, 0, 16, 10000)

    def test_physical_scales_must_be_positive(self):
        for P, XR in ((0, 1), (1, 0), (-1, 1)):
            with self.assertRaises(ValueError):
                physical_to_rows([Jet2(F(0))]*5, 0, P, XR)

    def test_nonfinite_scales_rejected(self):
        for P in (mp.inf, mp.nan):
            with self.assertRaises(ValueError):
                physical_to_rows([Jet2(F(0))]*5, 0, P, 10000)

    def test_angular_domain_enforced(self):
        for eta in (F(3, 2), mp.inf, mp.nan):
            with self.assertRaises(ValueError):
                physical_to_rows([Jet2(F(0))]*5, eta, 16, 10000)


class ReviewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.normalization = normalization_review()

    def test_all_rows_against_direct_differentiation(self):
        self.assertTrue(self.normalization['checks']['all_five_rows_keep_second_derivatives'])

    def test_normalization_curvature_cannot_be_frozen(self):
        self.assertTrue(self.normalization['checks']['frozen_normalization_failure_detected'])

    def test_C1_factor_is_not_a_C2_factor(self):
        self.assertTrue(self.normalization['checks']['old_factor_12_rejected'])

    def test_energy_cancellation(self):
        self.assertTrue(all(cancellation_review()['checks'].values()))

    def test_exponential_curvature(self):
        self.assertTrue(all(exponential_review()['checks'].values()))

    def test_pressure_axis_integrability(self):
        self.assertTrue(all(axis_integrability_review()['checks'].values()))


class BoundTests(unittest.TestCase):
    def test_actual_core_and_continuation_at_both_precisions(self):
        for digits in (80, 110):
            r = core_transfer(digits)
            self.assertTrue(all(r['checks'].values()))
            self.assertFalse(r['raw_C_or_width_materialized'])

    def test_all_five_normalized_rows(self):
        for digits in (80, 110):
            r = normalization_constants(digits)
            self.assertEqual(len(r['row_coefficients']), 5)
            self.assertEqual(r['inverse_f_squared_C2'], 20)
            self.assertTrue(all(r['checks'].values()))

    def test_actual_matching_C2_contraction(self):
        for digits in (80, 110):
            self.assertTrue(all(matching_transfer(digits)['checks'].values()))

    def test_raw_second_derivative_keeps_factorial(self):
        r = matching_transfer()
        self.assertGreater(lower(r['target_raw_second_derivative_bound']),
                           upper(r['target_C2_bound']))
        self.assertGreater(lower(r['coefficient_raw_second_from_norm']),
                           upper(r['coefficient_C2_bound']))

    def test_previous_C1_counterexample_is_still_retained(self):
        rows = previous.c1_counterexamples()
        self.assertTrue(all(r['C1_upper'] <= F(1, 10**16) for r in rows))
        self.assertGreater(rows[-1]['second_derivative_sup'], 1)

    def test_inherited_sources_unchanged(self):
        self.assertTrue(all(inherited_provenance().values()))


class ScopeTests(unittest.TestCase):
    def record(self):
        return deepcopy(dict(status=STATUS, flags=FLAGS, resolved_obligations=RESOLVED,
                             remaining_obligations=REMAINING, gates={'valid': True}))

    def test_only_actual_C2_transfer_closed(self):
        validate(self.record())
        self.assertEqual(RESOLVED, ['actual_target_C2_transfer'])
        self.assertEqual(len(REMAINING), 3)

    def test_frequency_promotion_rejected(self):
        r = self.record()
        r['flags']['actual_joined_profile_frequency_selected'] = True
        with self.assertRaises(ValueError):
            validate(r)

    def test_pressure_obligation_cannot_disappear(self):
        r = self.record()
        del r['remaining_obligations']['actual_modulation_error_constants']
        with self.assertRaises(ValueError):
            validate(r)

    def test_failed_estimate_rejected(self):
        r = self.record()
        r['gates']['valid'] = False
        with self.assertRaises(ValueError):
            validate(r)

    def test_changed_norm_or_tolerance_rejected(self):
        original = json.loads((HERE/'protocol.json').read_text())
        validate_protocol(original)
        for key, value in (('entry_epsilon', '1e-12'), ('angular_norm', 'C1')):
            changed = dict(original, **{key: value})
            with self.assertRaises(ValueError):
                validate_protocol(changed)


if __name__ == '__main__':
    unittest.main()
