#!/usr/bin/env python3
"""Higher derivatives, omitted-term controls and precise scientific scope."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp

# Import the current audit before legacy diagnostic imports modify sys.path.
from audit import (FLAGS, STATUS, RESOLVED, REMAINING, expected_protocol,
                   provenance, validate, compare_records, main)
from higher_bounds import (HERE, INPUT_POWER, transfer_C4, cutoff_certificate,
                           compact_exponents, scale_checks, inverse_exp,
                           loop_ledger, final_ledger, fourth_implicit_terms)
from higher_checks import (variance_review, root_review, phase_review, modulation_review,
                           fourth_implicit_review, repair_C3_review, moment_review,
                           derivative_gap_review, reference_review)
from higher_jets import regular_variance_factor, fixed_frequency_derivative


class HigherDerivativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.variance = variance_review()
        cls.root = root_review()
        cls.phase = phase_review()
        cls.modulation = modulation_review()
        cls.repair = repair_C3_review()
        cls.moments = moment_review()

    def test_regular_variance_positive_series_and_zero_limits(self):
        self.assertTrue(all(self.variance['checks'].values()))
        self.assertEqual(regular_variance_factor(0), mp.mpf('.5'))
        self.assertEqual(self.variance['variance_mu_at_zero'], '0.0')

    def test_implicit_third_root_includes_moving_pressure_and_zero_target(self):
        self.assertTrue(all(self.root['checks'].values()))
        self.assertEqual(self.root['samples'][1]['target_constant'], '0')

    def test_third_inverse_phase_and_pulled_back_primitives(self):
        self.assertTrue(all(self.phase['checks'].values()))
        self.assertEqual(len(self.phase['samples']), 6)
        self.assertTrue(all(mp.mpf(row['frozen_phase_gap']) > mp.mpf('1e-8') for row in self.phase['samples']))

    def test_fixed_frequency_chain_keeps_growing_radial_term(self):
        self.assertTrue(all(self.modulation['checks'].values()))
        self.assertEqual(self.modulation['chain_sample_count'], 15)
        self.assertGreater(mp.mpf(self.modulation['missing_N_phase_second_gap']), 1)

    def test_exponential_third_derivative_keeps_cube_and_cross_terms(self):
        self.assertLess(mp.mpf(self.modulation['exponential_error']), mp.mpf('1e-55'))
        self.assertGreater(mp.mpf(self.modulation['missing_exponential_cube_gap']), mp.mpf('1e-4'))

    def test_final_nonlinear_C3_repair_and_normalization(self):
        self.assertTrue(all(self.repair['checks'].values()))
        self.assertTrue(self.repair['quadrature_is_not_a_continuum_certificate'])

    def test_original_fourth_moments_use_rational_integrals(self):
        self.assertTrue(all(self.moments['checks'].values()))
        self.assertEqual(len(self.moments['samples']), 5)
        self.assertTrue(self.moments['axis_pressure_integrand_uses_F_squared'])

    def test_all_eight_fourth_implicit_terms_are_retained(self):
        review = fourth_implicit_review()
        self.assertTrue(all(review['checks'].values()))
        self.assertEqual(fourth_implicit_terms()['coefficients'], [8, 6, 8, 24, 12, 12, 8, 1])

    def test_C3_cannot_silently_supply_the_C4_input(self):
        review = derivative_gap_review()
        self.assertTrue(all(review['checks'].values()))
        self.assertEqual(review['rows'][-1]['raw_fourth_peak'], 512)

    def test_reference_source_integral_matches_original_moments_through_C4(self):
        review = reference_review()
        self.assertTrue(all(review['checks'].values()))
        self.assertGreater(mp.mpf(review['frozen_axis_gradient_gap']), mp.mpf('1e-4'))


class ContinuumBoundTests(unittest.TestCase):
    def test_reference_C4_is_derived_and_compared_at_two_precisions(self):
        for digits in (80, 110):
            record = transfer_C4(digits)
            self.assertTrue(all(record['checks'].values()))
            self.assertTrue(record['reference_fourth_order_derived_from_defining_equations'])
            self.assertEqual(record['exact_C4_evaluation'], '14624/6075')
            self.assertEqual(record['normalized_energy_row_coefficient'], 8125)
            self.assertEqual(record['f_inverse_square_C4'], 25)

    def test_cutoff_edges_and_inverse_widths_have_a_uniform_fifth_bound(self):
        for digits in (80, 110):
            record = cutoff_certificate(digits)
            self.assertEqual(record['maximum_step_order'], 5)
            self.assertTrue(all(record['checks'].values()))

    def test_mixed_inverse_jet_cost_and_higher_input_padding(self):
        record = compact_exponents()
        self.assertEqual(inverse_exp(2048, 4, 4), 8213)
        self.assertTrue(all(record['checks'].values()))
        self.assertEqual(record['higher_input_A_power'], 256)
        self.assertTrue(record['original_value_bounds_and_parameters_unchanged'])

    def test_same_N_and_H_absorb_higher_derivatives(self):
        for digits in (80, 110):
            record = scale_checks(digits)
            self.assertTrue(all(record['checks'].values()))
            self.assertEqual(record['fixed_frequency'], '1+floor(H^32)')
            self.assertFalse(record['large_parameters_or_frequency_materialized'])
        self.assertTrue(all(loop_ledger()['checks'].values()))

    def test_explicit_K_still_gives_singular_background_residual_bounds(self):
        record = final_ledger()
        self.assertEqual(record['ordinary_jet_envelope'], 'K=H^64')
        self.assertFalse(record['residual_is_proved_to_decay'])
        self.assertFalse(record['full_wave_tensor_or_global_physical_residual_included'])
        self.assertTrue(all(record['checks'].values()))

    def test_invalid_inverse_orders_and_frequencies_are_rejected(self):
        for values in ((1, -1, 2), (1, 2, 1.5)):
            with self.assertRaises(ValueError):
                inverse_exp(*values)
        with self.assertRaises(ValueError):
            fixed_frequency_derivative(lambda y, eta, s: y, 0, 0, 0, 1, 0)


class ProvenanceAndScopeTests(unittest.TestCase):
    def record(self):
        return deepcopy(dict(status=STATUS, flags=FLAGS, resolved_obligations=RESOLVED,
                             remaining_obligations=REMAINING, gates={'valid': True}))

    def test_protocol_and_all_preceding_manuscript_hashes(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())
        old = provenance()
        self.assertGreater(len(old), 261)
        self.assertTrue(all(old.values()))

    def test_global_solution_and_blowup_promotions_are_rejected(self):
        validate(self.record())
        for key in ('actual_global_first_coefficient_extension_verified', 'full_physical_residual_budget_verified',
                    'full_physical_stress_wave_construction_verified', 'physical_background_residual_decay_verified',
                    'full_PDE_corrections_verified', 'smooth_force_verified', 'blowup_verified'):
            record = self.record(); record['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(record)

    def test_omitted_obligations_and_failed_controls_are_rejected(self):
        for field, value in (('remaining_obligations', []), ('gates', {'failed': False})):
            record = self.record(); record[field] = value
            with self.assertRaises(ValueError):
                validate(record)

    def test_review_roundoff_does_not_relax_metadata_or_large_errors(self):
        compare_records({'error': '1e-75'}, {'error': '2e-75'}, 'root.reviews.test')
        for location in ('root.protocol', 'root.reviews.test'):
            with self.assertRaises(ValueError):
                compare_records({'error': '1e-75'}, {'error': '1e-40'}, location)
        with self.assertRaises(ValueError):
            compare_records({'N': 'H^32'}, {'N': 'H^33'}, 'root.protocol')

    def test_existing_evidence_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'keep.json'; path.write_text('keep')
            with self.assertRaises(SystemExit):
                main(['--output', str(path)])
            self.assertEqual(path.read_text(), 'keep')


if __name__ == '__main__':
    unittest.main()
