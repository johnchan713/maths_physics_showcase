#!/usr/bin/env python3
"""Physical identities, missing-term controls and limits on proof promotion."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from residual import (HERE, derivative_loss, evaluate_series, minimum_order_for_decay,
                      momentum, physical_point, physical_series, series_coefficients)
from bounds import ledger
from review import (cutoff_review, derivative_gap_review, exponent_review, fixtures,
                    omission_review, physical_review, to_cartesian)
from audit import (FLAGS, REMAINING, RESOLVED, STATUS, compare_records,
                   expected_protocol, main, provenance, validate)


class PhysicalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.physical = physical_review()
        cls.omissions = omission_review()

    def test_complete_base_residual_against_cartesian_jets(self):
        self.assertLess(self.physical['maxima']['base'], 2e-11)
        self.assertEqual(self.physical['physical_sample_count'], 324)
        self.assertEqual(self.physical['maxima']['old_ansatz'], 0.)

    def test_positive_order_all_products_and_divergence(self):
        self.assertLess(self.physical['maxima']['finite_correction'], 2e-11)
        self.assertLess(self.physical['maxima']['divergence'], 2e-11)
        self.assertTrue(self.omissions['checks']['missing_positive_order_lambda'])

    def test_stress_pair_and_full_cross_tensor_are_distinct(self):
        self.assertLess(self.physical['maxima']['leading_stress'], 2e-11)
        self.assertLess(self.physical['maxima']['full_symmetric_stress'], 2e-11)
        self.assertGreater(self.omissions['negative_control_gaps']['missing_symmetric_rz_axial_divergence'], 1e-4)

    def test_physical_scalar_stencil_refines(self):
        errors = self.omissions['finite_difference_errors']
        self.assertLess(errors[-1], 1e-6)
        self.assertGreater(errors[0], 4*errors[-1])

    def test_all_six_missing_term_controls_fail(self):
        self.assertEqual(len(self.omissions['negative_control_gaps']), 6)
        self.assertTrue(all(v > 1e-4 for v in self.omissions['negative_control_gaps'].values()))

    def test_three_coefficient_orders_retain_all_convolutions(self):
        # Additional order not used in the reported two-coefficient fixtures.
        base, correction = fixtures()
        third = deepcopy(correction)
        point = physical_point(.17, -.4, .8, -.7, .11)
        u, p = physical_series([base, correction, third], point, .11)
        actual, scales = momentum(u, p, .9)
        coefficients = series_coefficients([base, correction, third], .8, -.4, .11, .9)
        expected = to_cartesian(evaluate_series(coefficients, .17, .8, .11), -.7)
        self.assertLess(float(np.max(np.abs(actual-expected)/scales)), 2e-11)
        self.assertEqual(len(coefficients['radial']), 6)
        self.assertEqual(len(coefficients['theta']), 5)


class BudgetTests(unittest.TestCase):
    def test_C2_family_has_unbounded_radial_axial_diffusion(self):
        review = derivative_gap_review()
        self.assertTrue(all(review['checks'].values()))
        self.assertEqual(review['C2_counterexample'][0]['radial_axial_diffusion_coefficient'], '49/8')

    def test_cutoffs_act_on_streamfunctions_to_preserve_divergence(self):
        review = cutoff_review()
        self.assertTrue(all(review['checks'].values()))
        self.assertFalse(review['diagnostic_cutoff_is_compact'])

    def test_physical_derivative_losses_keep_their_distinct_exponents(self):
        self.assertEqual(derivative_loss(F(1, 200), 2, 1, 1), F(499, 200))
        self.assertEqual(minimum_order_for_decay(F(1, 200)), 151)
        self.assertEqual(minimum_order_for_decay(F(1, 200), time=2), 351)
        self.assertTrue(all(exponent_review()['checks'].values()))
        for counts in ((-1, 0, 0), (0, 1.5, 0)):
            with self.assertRaises(ValueError):
                derivative_loss(F(1, 200), *counts)

    def test_annular_constants_are_conditional_and_outward(self):
        record = ledger()
        self.assertEqual(record['axial_second_polynomial'], [22, 18, 4])
        self.assertTrue(all(record['checks'].values()))
        self.assertFalse(record['actual_K_evaluated'])

    def test_frozen_protocol_and_preceding_sources(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())
        self.assertTrue(all(provenance().values()))


class ScopeTests(unittest.TestCase):
    def record(self):
        return deepcopy(dict(status=STATUS, flags=FLAGS, resolved_obligations=RESOLVED,
                             remaining_obligations=REMAINING, gates={'valid': True}))

    def test_actual_C3_and_positive_order_solution_remain_open(self):
        validate(self.record())
        self.assertTrue(FLAGS['inherited_C2_does_not_bound_full_residual'])
        for key in ('actual_post_modulation_C3_bound_verified', 'first_correction_for_actual_profile_solved'):
            record = self.record()
            record['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(record)

    def test_wave_PDE_force_and_blowup_promotions_rejected(self):
        for key in ('full_physical_residual_budget_verified', 'full_physical_stress_wave_construction_verified',
                    'full_PDE_corrections_verified', 'smooth_force_verified', 'blowup_verified'):
            record = self.record()
            record['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(record)

    def test_missing_obligations_or_failed_controls_rejected(self):
        record = self.record()
        record['remaining_obligations'] = []
        with self.assertRaises(ValueError):
            validate(record)
        record = self.record()
        record['gates']['valid'] = False
        with self.assertRaises(ValueError):
            validate(record)

    def test_roundoff_tolerance_does_not_allow_status_changes(self):
        compare_records({'a': 1e-14, 'b': False}, {'a': 2e-14, 'b': False})
        for record in ({'a': .1, 'b': False}, {'a': 1e-14, 'b': True}):
            with self.assertRaises(ValueError):
                compare_records({'a': 1e-14, 'b': False}, record)

    def test_existing_outputs_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'keep.json'
            path.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output', str(path)])
            self.assertEqual(path.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
