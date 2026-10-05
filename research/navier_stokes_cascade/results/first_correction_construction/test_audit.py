#!/usr/bin/env python3
"""Test distinct correction identities, inherited bounds and scope limits."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
import numpy as np

from transfer import HERE, angular_transfer, inner_operator_certificate, norm_ledger, third_implicit_terms
from majorant import allowed_word, log_term, log_sum_bound
from system import matrices
from series import build
from checks import (matrix_review, sparse_review, recurrence_review,
                    first_moment_review, implicit_review, fixtures)
from audit import (FLAGS, STATUS, RESOLVED, REMAINING, expected_protocol, provenance,
                   validate, compare_records, main, protected_growth_ledger)


class CorrectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.matrix = matrix_review()
        cls.recurrence = recurrence_review()
        cls.sparse = sparse_review()
        cls.moments = first_moment_review()

    def test_system_against_independent_physical_coefficient_assembly(self):
        self.assertEqual(self.matrix['sample_count'], 54)
        self.assertLess(self.matrix['maximum_relative_error'], 2e-12)

    def test_pressure_lambda_and_angular_omissions_are_detected(self):
        self.assertTrue(all(value > 1e-5 for value in self.matrix['negative_controls'].values()))

    def test_sparse_matrix_cancels_also_when_derivative_hits_coefficients(self):
        base, _ = fixtures()
        left = matrices(base, .2, .3, .005)[1]
        diff = (matrices(base, .4, .3+1e-5, .005)[1]-matrices(base, .4, .3-1e-5, .005)[1])/2e-5
        self.assertEqual(float(np.max(np.abs(left@np.diag(np.arange(1., 7.))@diff))), 0.)
        self.assertTrue(self.matrix['checks']['radial_parity'])

    def test_word_count_controls_angular_loss_without_small_radius(self):
        self.assertTrue(all(self.sparse['checks'].values()))
        self.assertFalse(allowed_word((1, 1)))
        self.assertTrue(allowed_word((1, 0, 1, 0, 1)))
        self.assertGreater(100*float('.1'), 1)
        self.assertTrue(self.sparse['finite_sum_checks'][-1]['below_bound'])

    def test_first_radial_slopes_and_degree_refinement(self):
        self.assertTrue(all(self.recurrence['checks'].values()))
        self.assertEqual([row['degree'] for row in self.recurrence['rows']], [3, 5, 7])
        self.assertLess(max(self.recurrence['rows'][-1]['errors']), 2e-11)

    def test_axis_traces_are_zero_and_angular_budget_is_retained(self):
        base, _ = fixtures()
        correction, budget = build(base, .005, 1., 3)
        for eta in (-.5, 0., .5):
            for field in (correction.F, correction.U, correction.Pi, correction.average_U):
                self.assertEqual(field(0., eta), 0.)
        self.assertEqual(budget['retained_angular_orders'], [23, 22, 21])

    def test_positive_order_moments_are_linear_for_arbitrary_finite_target(self):
        self.assertTrue(all(self.moments['checks'].values()))
        self.assertEqual(self.moments['positive_order_patch'], 'I_3')
        self.assertEqual(self.moments['rows'][-1]['lambda_value'], '1e-50')
        self.assertFalse(self.moments['actual_annular_support_hypotheses_fully_verified'])

    def test_third_implicit_derivative_retains_all_five_product_terms(self):
        self.assertTrue(all(implicit_review()['checks'].values()))
        self.assertEqual(third_implicit_terms()['coefficients'], [6, 6, 6, 6, 1])
        self.assertEqual(third_implicit_terms()['raw_derivative_factor'], 6)


class ActualBoundTests(unittest.TestCase):
    def test_unmodulated_C3_transfer_at_two_interval_precisions(self):
        for digits in (80, 110):
            record = angular_transfer(digits)
            self.assertTrue(all(record['checks'].values()))
            self.assertEqual(record['exact_evaluation_constant'], '544/243')
            self.assertFalse(record['post_modulation_C3_bound_supplied'])
            self.assertFalse(record['global_third_order_compact_envelope_supplied'])

    def test_C3_normalization_and_missing_C4_input(self):
        ledger = norm_ledger()
        self.assertEqual(ledger['f_inverse_square_C3'], Fraction(24))
        self.assertEqual(ledger['energy_row_coefficient'], 7800)
        self.assertEqual(ledger['required_original_moment_order_for_that_loop'], 4)

    def test_actual_inner_bound_has_strict_core_domain_and_no_L_pole(self):
        for digits in (80, 110):
            record = inner_operator_certificate(digits)
            self.assertTrue(all(record['checks'].values()))
            self.assertEqual(record['actual_operator_envelope'], 'C^16')
            self.assertFalse(record['actual_inner_interval_numerically_resolved'])

    def test_majorant_refuses_unproved_zero_bounds_or_missing_geometry(self):
        for C, a, delta in ((0, 1, .1), (1, 0, .1), (1, 1, 0), (1, 1, 2)):
            with self.assertRaises(ValueError):
                log_term(2, C, a, delta)
            with self.assertRaises(ValueError):
                log_sum_bound(C, a, delta)

    def test_protected_partial_growth_does_not_claim_the_complete_field(self):
        ledger = protected_growth_ledger()
        self.assertEqual(ledger['point'], 'X=1/Lambda, eta=0')
        self.assertFalse(ledger['threshold_materialized'])
        self.assertFalse(ledger['complete_corrected_field_bound_supplied'])


class ProvenanceScopeTests(unittest.TestCase):
    def record(self):
        return deepcopy(dict(status=STATUS, flags=FLAGS, resolved_obligations=RESOLVED,
                             remaining_obligations=REMAINING, gates={'valid': True}))

    def test_protocol_and_frozen_manuscript_sources(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())
        self.assertTrue(all(provenance().values()))

    def test_global_extension_waves_force_and_blowup_cannot_be_promoted(self):
        validate(self.record())
        for key in ('actual_post_modulation_C3_bound_verified', 'actual_global_first_coefficient_extension_verified',
                    'full_physical_stress_wave_construction_verified', 'full_PDE_corrections_verified',
                    'smooth_force_verified', 'complete_field_energy_verified', 'blowup_verified'):
            record = self.record(); record['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(record)

    def test_removing_a_gap_or_ignoring_failed_controls_is_rejected(self):
        for field, value in (('remaining_obligations', []), ('gates', {'failed': False})):
            record = self.record(); record[field] = value
            with self.assertRaises(ValueError):
                validate(record)

    def test_frozen_numeric_comparison_accepts_roundoff_only(self):
        compare_records({'x': 1e-14, 'closed': False}, {'x': 2e-14, 'closed': False})
        for value in ({'x': .1, 'closed': False}, {'x': 1e-14, 'closed': True}):
            with self.assertRaises(ValueError):
                compare_records({'x': 1e-14, 'closed': False}, value)

    def test_existing_evidence_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'keep.json'; path.write_text('keep')
            with self.assertRaises(SystemExit):
                main(['--output', str(path)])
            self.assertEqual(path.read_text(), 'keep')


if __name__ == '__main__':
    unittest.main()
