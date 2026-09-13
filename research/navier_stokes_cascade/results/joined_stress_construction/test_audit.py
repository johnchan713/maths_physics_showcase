#!/usr/bin/env python3
"""Tests target the new derivative and representation risks, not blow-up."""
from copy import deepcopy
from fractions import Fraction
import unittest
import mpmath as mp
from construction import (CANDIDATE, OPEN_OBLIGATIONS, geometry, derivative_constants,
                          curvature_bound, c1_counterexamples, finite_integer_comparison,
                          promotion_review, lower, upper)
from review import curvature_identity, loop_derivative_controls, rounded_endpoint_control
from audit import STATUS, validate, inherited_hashes


class GeometryTests(unittest.TestCase):
    def test_actual_geometry_at_both_precisions(self):
        for digits in (80, 110):
            g = geometry(digits)
            self.assertTrue(all(g['checks'].values()))
            self.assertFalse(g['endpoint_offsets_materialized'])

    def test_rounded_endpoints_cannot_certify_separation(self):
        self.assertTrue(rounded_endpoint_control()['rounded_offsets_collapse'])

    def test_frozen_provenance(self):
        self.assertTrue(all(inherited_hashes().values()))


class DerivativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.identity = curvature_identity()

    def test_original_map_second_derivative(self):
        self.assertTrue(self.identity['checks']['second_derivative_identity'])

    def test_mixed_term_cannot_be_dropped(self):
        self.assertTrue(self.identity['checks']['missing_cross_detected'])

    def test_explicit_angular_coefficient_cannot_be_frozen(self):
        self.assertTrue(self.identity['checks']['frozen_zeta_detected'])

    def test_second_derivative_and_cutoff_constants(self):
        for digits in (80, 110):
            self.assertTrue(all(derivative_constants(digits)['checks'].values()))

    def test_conditional_curvature_transfer_retains_Z2(self):
        self.assertGreater(lower(curvature_bound('1e4')), 1000*10000)
        self.assertLess(upper(curvature_bound('1e4')), 1001*10000)
        self.assertGreater(lower(curvature_bound(0)), 0)

    def test_invalid_curvature_bound_rejected(self):
        for value in ('-1', 'nan', 'inf'):
            with self.assertRaises(ValueError):
                curvature_bound(value)

    def test_C1_cannot_be_promoted_to_C2(self):
        rows = c1_counterexamples()
        self.assertTrue(all(r['C1_upper'] <= Fraction(1, 10**16) for r in rows))
        self.assertTrue(all(a['value_sup'] > b['value_sup'] and
                            a['second_derivative_sup'] < b['second_derivative_sup']
                            for a, b in zip(rows, rows[1:])))
        self.assertEqual(rows[-1]['second_derivative_sup'], 50000000)

    def test_loop_root_derivative_including_zero_pressure(self):
        self.assertTrue(loop_derivative_controls()['check'])

    def test_cap_coefficient_inequalities(self):
        self.assertTrue(all(finite_integer_comparison().values()))


class PromotionTests(unittest.TestCase):
    def record(self):
        return deepcopy(dict(status=STATUS, promotion=promotion_review(),
                             proposed_construction=CANDIDATE,
                             open_obligations=OPEN_OBLIGATIONS,
                             gates={'valid_control': True}))

    def test_review_accepts_record_of_unclosed_work(self):
        validate(self.record())

    def test_failed_check_rejected(self):
        r = self.record()
        r['gates']['valid_control'] = False
        with self.assertRaises(ValueError):
            validate(r)

    def test_unproved_frequency_promotion_rejected(self):
        r = self.record()
        r['promotion']['actual_joined_profile_frequency_selected'] = True
        with self.assertRaises(ValueError):
            validate(r)

    def test_missing_obligation_cannot_be_erased(self):
        r = self.record()
        del r['open_obligations']['actual_target_C2_transfer']
        with self.assertRaises(ValueError):
            validate(r)

    def test_proposal_cannot_be_relabelled_a_certificate(self):
        r = self.record()
        r['proposed_construction']['certified_for_actual_joined_profile'] = True
        with self.assertRaises(ValueError):
            validate(r)


if __name__ == '__main__':
    unittest.main()
