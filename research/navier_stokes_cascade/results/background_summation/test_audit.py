#!/usr/bin/env python3
"""General-order sources, actual curl product rules, summation and scope."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
import numpy as np

# Load this audit before inherited modules alter the import search path.
from audit import (STATUS, FLAGS, RESOLVED, REMAINING, expected_protocol,
                   provenance, validate, compare_records, main)
from hierarchy import HERE, matrices, lower_sources, induction_ledger, finite_residual_exponent
from summation import cutoff_log_scale, scalar_coefficients, scalar_residual_coefficients, loss_budget, proof_ledger
from checks import (fixtures, matrix_review, moment_review, cutoff_review, scale_review,
                    scalar_summation_review, exponent_review)


class BackgroundTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.matrices = matrix_review()
        cls.moments = moment_review()
        cls.cutoffs = cutoff_review()
        cls.scales = scale_review()
        cls.scalar = scalar_summation_review()
        cls.exponents = exponent_review()

    def test_general_system_matches_separate_complete_residual(self):
        self.assertEqual(self.matrices['sample_count'], 216)
        self.assertEqual(self.matrices['tested_orders'], [1, 2, 3, 5])
        self.assertTrue(self.matrices['checks']['independent_complete_residual_system'])

    def test_frozen_first_order_system_is_preserved(self):
        self.assertTrue(self.matrices['checks']['first_order_system_preserved'])

    def test_shifted_radial_axial_viscosity_enters_pressure_at_order_two(self):
        self.assertTrue(self.matrices['checks']['shifted_radial_viscosity_detected'])
        self.assertEqual(self.matrices['controls'][0]['order'], 2)
        for row in self.matrices['controls']:
            self.assertAlmostEqual(row['omitted_radial_axial_pressure_gap'], row['expected_half_axial_source'], places=12)

    def test_sources_are_finite_at_the_axis_for_every_tested_order(self):
        for n in (1, 2, 3, 5):
            self.assertTrue(all(np.isfinite(a).all() for a in matrices(fixtures(n), 0., .4, .005)))

    def test_derivative_matrix_and_its_sparse_product_are_order_independent(self):
        reference = matrices(fixtures(1), .4, .3, .005)[1]
        for n in (2, 3, 5):
            self.assertTrue(np.array_equal(matrices(fixtures(n), .4, .3, .005)[1], reference))
        delta = (matrices(fixtures(3), .5, .3+1e-5, .005)[1]-matrices(fixtures(3), .5, .3-1e-5, .005)[1])/2e-5
        self.assertEqual(float(np.max(abs(reference@np.diag(np.arange(1., 7.))@delta))), 0.)

    def test_all_order_moment_identities_retain_nonzero_lower_viscosity(self):
        self.assertEqual(self.moments['sample_count'], 48)
        self.assertTrue(all(self.moments['checks'].values()))

    def test_cutoff_curl_is_divergence_free_at_axis_and_transition(self):
        self.assertEqual(self.cutoffs['sample_count'], 96)
        self.assertTrue(self.cutoffs['includes_axis_and_both_plateaus'])
        self.assertTrue(self.cutoffs['checks']['exact_curl_cutoff_is_divergence_free'])

    def test_independently_differentiated_potential_reproduces_cut_velocity(self):
        self.assertTrue(self.cutoffs['checks']['separate_potential_curl_agrees'])

    def test_missing_or_wrong_sign_commutators_are_rejected(self):
        self.assertTrue(self.cutoffs['checks']['direct_velocity_cutoff_rejected'])
        self.assertTrue(self.cutoffs['checks']['wrong_commutator_sign_rejected'])

    def test_log_scales_control_extreme_requirements_at_two_precisions(self):
        self.assertEqual([x['digits'] for x in self.scales['records']], [80, 110])
        self.assertTrue(all(self.scales['checks'].values()))
        self.assertFalse(self.scales['actual_NS_coefficient_bounds_or_scales_evaluated'])

    def test_scale_selector_rejects_missing_or_invalid_proved_bounds(self):
        for index, gain, constraints in ((0, 1, [(0, 0)]), (1, 0, [(0, 0)]), (1, 1, []),
                                         (1, 1, [(-1, 0)]), (1, 1, [(0, -1)]), (1, 'inf', [(0, 0)])):
            with self.assertRaises(ValueError):
                cutoff_log_scale(index, gain, constraints)

    def test_scalar_formal_series_has_factorial_growth_and_exact_recurrence(self):
        self.assertEqual(scalar_coefficients(4), [0, 1, 2, 8, 44])
        self.assertTrue(self.scalar['checks']['formal_recurrence_cancels_all_retained_orders'])
        self.assertTrue(self.scalar['checks']['factorial_growth'])
        self.assertTrue(self.scalar['formal_series_is_divergent_at_every_positive_q'])

    def test_scalar_cutoff_sum_is_locally_finite_with_improving_residual(self):
        self.assertTrue(all(self.scalar['checks'].values()))
        self.assertEqual([r['active_index'] for r in self.scalar['rows']], [4, 8, 12, 16])
        self.assertTrue(self.scalar['scalar_ODE_only'])

    def test_finite_exponents_improve_without_order_dependent_derivative_loss(self):
        self.assertTrue(all(self.exponents['checks'].values()))
        self.assertEqual(finite_residual_exponent(Fraction(1, 200), 1), Fraction(-149, 100))
        self.assertEqual(finite_residual_exponent(Fraction(1, 200), 150), 0)
        self.assertGreater(finite_residual_exponent(Fraction(1, 200), 151), 0)

    def test_induction_has_common_support_but_no_uniform_analytic_strip_claim(self):
        ledger = induction_ledger()
        self.assertTrue(ledger['angular_strip_may_shrink_with_order'])
        self.assertTrue(ledger['no_uniform_coefficient_bound_in_order_asserted'])
        self.assertTrue(ledger['higher_stress_over_zeta_has_all_finite_derivative_bounds'])

    def test_background_flatness_and_growth_exclude_physical_wave_field(self):
        ledger = proof_ledger()
        self.assertTrue(ledger['augmented_background_residual_flat_in_all_fixed_physical_derivatives'])
        self.assertFalse(ledger['bare_background_momentum_residual_flat'])
        self.assertFalse(ledger['lower_bound_includes_physical_waves'])
        self.assertTrue(ledger['uniform_perturbed_stress_cone_not_inferred'])

    def test_invalid_orders_and_scalar_inputs_are_rejected(self):
        for args in ((0, 1, 0), (Fraction(1, 100), 1, 0), (Fraction(1, 200), -1, 0)):
            with self.assertRaises(ValueError):
                finite_residual_exponent(*args)
        with self.assertRaises(ValueError):
            loss_budget(Fraction(1, 200), -1)
        with self.assertRaises(ValueError):
            scalar_coefficients(0)
        with self.assertRaises(ValueError):
            scalar_residual_coefficients([0, 1], [0], [0])


class ProvenanceScopeTests(unittest.TestCase):
    def record(self):
        return deepcopy(dict(status=STATUS, flags=FLAGS, resolved_obligations=RESOLVED,
                             remaining_obligations=REMAINING, gates={'valid': True}))

    def test_protocol_and_all_prior_source_hashes(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())
        inherited = provenance()
        self.assertGreater(len(inherited), 293)
        self.assertTrue(all(inherited.values()))

    def test_complete_solution_and_bare_residual_promotions_are_rejected(self):
        validate(self.record())
        for key in ('bare_background_momentum_residual_flat', 'uniform_corrected_stress_cone_verified',
                    'full_physical_stress_wave_construction_verified', 'full_physical_residual_budget_verified',
                    'smooth_force_verified', 'complete_field_energy_verified', 'complete_field_blowup_lower_bound_verified', 'blowup_verified'):
            record = self.record(); record['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(record)

    def test_failed_controls_and_removed_obligations_are_rejected(self):
        for field, value in (('gates', {'failed': False}), ('remaining_obligations', [])):
            record = self.record(); record[field] = value
            with self.assertRaises(ValueError):
                validate(record)

    def test_frozen_rationals_and_scope_are_not_rounded_away(self):
        compare_records({'value': '1e-75'}, {'value': '2e-75'}, 'root.reviews.test')
        for a, b in (('1/200', '1/201'), ('1e-75', '1e-30')):
            with self.assertRaises(ValueError):
                compare_records({'value': a}, {'value': b}, 'root.reviews.test')
        with self.assertRaises(ValueError):
            compare_records({'N': 'H^32'}, {'N': 'H^33'}, 'root.protocol')

    def test_existing_evidence_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'keep.json'; path.write_text('keep')
            with self.assertRaises(SystemExit):
                main(['--output', str(path)])
            self.assertEqual(path.read_text(), 'keep')


if __name__ == '__main__':
    unittest.main()
