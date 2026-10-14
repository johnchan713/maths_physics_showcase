#!/usr/bin/env python3
"""Test the corrected geometry, global moments, retained source and scope."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp

# Import this audit before inherited numerical imports alter sys.path.
from audit import (STATUS, FLAGS, RESOLVED, REMAINING, expected_protocol,
                   provenance, validate, compare_records, main)
from extension import (HERE, normalized_target, solve_patch, geometry_ledger,
                       construction_ledger, residual_gain_ledger, Repair)
from checks import (kinematic_review, global_matching_review, memory_source_review,
                    moment_identity_review, remainder_review)
from fixtures import MemoryFixture, Rule


class ExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.global_moments = global_matching_review()
        cls.source = memory_source_review()
        cls.moment_identities = moment_identity_review()
        cls.remainder = remainder_review()

    def test_cutoff_reaches_outside_core_before_modulation(self):
        record = geometry_ledger()
        self.assertTrue(all(record['checks'].values()))
        self.assertEqual(record['inner_offsets_in_units_tc'],
                         list(map(Fraction, ('1/16', '1/8', '3/16', '1/4', '1/2'))))

    def test_I3_and_mean_patch_are_separated_and_pressure_edge_is_later(self):
        record = geometry_ledger()
        self.assertEqual(record['positive_patch_stage_offsets'], [-14, -9])
        self.assertEqual(record['mean_patch_stage_offsets'], [-8, -3])
        self.assertIn('Xv', record['global_pressure_support_edge'])

    def test_collar_bound_does_not_reuse_core_only_C16(self):
        record = construction_ledger()
        self.assertFalse(record['collar_bound_is_C16'])
        self.assertFalse(record['collar_bound_uniform_in_large_parameters'])
        self.assertTrue(record['same_inner_coefficient_by_analytic_uniqueness'])

    def test_divergence_retains_D_plus_lambda_exactly(self):
        self.assertTrue(all(kinematic_review()['checks'].values()))

    def test_nonzero_radial_memory_is_checked_independently(self):
        self.assertEqual(self.source['sample_count'], 15)
        self.assertTrue(all(self.source['checks'].values()))

    def test_twelve_global_manufactured_repairs_close_five_moments(self):
        self.assertEqual(len(self.global_moments['rows']), 12)
        self.assertTrue(self.global_moments['checks']['five_global_moments_close'])
        self.assertFalse(self.global_moments['actual_huge_frequency_profile_evaluated'])

    def test_omitted_source_changes_pressure_and_fifth_moments(self):
        self.assertTrue(self.global_moments['checks']['missing_pressure_source_detected_off_midplane'])
        self.assertTrue(self.global_moments['checks']['missing_fifth_source_detected'])

    def test_midplane_total_cancellation_does_not_erase_local_pressure(self):
        self.assertTrue(self.global_moments['checks']['midplane_total_pressure_source_cancels'])
        self.assertTrue(self.global_moments['checks']['pressure_need_not_vanish_on_mean_patch'])

    def test_patch_source_is_negative_even_when_U0_is_zero(self):
        fixture = MemoryFixture()
        for eta in map(mp.mpf, ('-.6', '0', '.6')):
            X = fixture.radius*mp.exp(2)
            self.assertEqual(fixture.memory_jets(X), [1, 0, 0, 0])
            self.assertAlmostEqual(fixture.omega(X, eta), -1/(2*X), places=65)

    def test_linear_inverse_accepts_large_targets_without_a_small_root(self):
        repair = Repair('.001', order=16, panels=2)
        moments = list(map(mp.mpf, ('1e8', '-3e7', '2e5', '7e8', '-4e9')))
        c = solve_patch(moments, 3, mp.mpf('.2'), repair)
        z = normalized_target(moments, 3, mp.mpf('.2'), repair.lam)
        self.assertLess(max(map(abs, repair.B*c+z)), mp.mpf('1e-55'))
        self.assertGreater(max(map(abs, c)), 1)

    def test_invalid_radius_amplitude_slope_and_nonfinite_data_rejected(self):
        for radius, amplitude, lam in ((0, 1, '.001'), (1, 0, '.001'), (1, 1, 0), ('inf', 1, '.001')):
            with self.assertRaises(ValueError):
                normalized_target([0]*5, radius, amplitude, lam)
        with self.assertRaises(ValueError):
            normalized_target([0, 0, 0, 0, 'nan'], 1, 1, '.001')
        with self.assertRaises(ValueError):
            normalized_target([0]*4, 1, 1, '.001')

    def test_conservative_moment_identities_use_nonzero_discrepancies(self):
        self.assertEqual(self.moment_identities['sample_count'], 18)
        self.assertTrue(all(self.moment_identities['checks'].values()))

    def test_full_quadratic_and_coefficient_diffusion_remainder(self):
        self.assertEqual(self.remainder['sample_count'], 18)
        self.assertTrue(all(self.remainder['checks'].values()))
        self.assertTrue(self.remainder['first_two_coefficients_not_asserted_zero_for_these_fixtures'])

    def test_relative_gain_is_not_absolute_decay_or_full_wave_error(self):
        record = residual_gain_ledger()
        self.assertEqual(record['relative_gain'], Fraction(1, 100))
        self.assertEqual(record['radial_physical_exponent'], Fraction(-149, 100))
        self.assertFalse(record['absolute_error_proved_to_decay'])
        self.assertFalse(record['error_is_full_wave_or_localized_residual'])
        self.assertEqual(residual_gain_ledger(derivative_order=3)['sufficient_all_component_exponent'], Fraction(-449, 100))

    def test_invalid_derivative_budget_is_rejected(self):
        for h, order in ((0, 0), (Fraction(1, 100), 0), (Fraction(1, 200), -1), (Fraction(1, 200), .5)):
            with self.assertRaises(ValueError):
                residual_gain_ledger(h, order)


class ProvenanceScopeTests(unittest.TestCase):
    def record(self):
        return deepcopy(dict(status=STATUS, flags=FLAGS, resolved_obligations=RESOLVED,
                             remaining_obligations=REMAINING, gates={'valid': True}))

    def test_protocol_and_inherited_bytes_preserved(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())
        inherited = provenance()
        self.assertGreater(len(inherited), 278)
        self.assertTrue(all(inherited.values()))

    def test_conditional_extension_cannot_promote_full_solution(self):
        validate(self.record())
        for key in ('full_physical_stress_wave_construction_verified', 'absolute_physical_residual_decay_verified',
                    'all_order_background_induction_verified', 'full_PDE_corrections_verified',
                    'smooth_force_verified', 'complete_field_energy_verified', 'blowup_verified'):
            record = self.record(); record['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(record)

    def test_failed_controls_or_removed_obligations_rejected(self):
        for field, value in (('gates', {'failed': False}), ('remaining_obligations', [])):
            record = self.record(); record[field] = value
            with self.assertRaises(ValueError):
                validate(record)

    def test_frozen_comparison_retains_tight_scope_and_numeric_checks(self):
        compare_records({'error': '1e-75'}, {'error': '2e-75'}, 'root.reviews.test')
        for where in ('root.protocol', 'root.reviews.test'):
            with self.assertRaises(ValueError):
                compare_records({'error': '1e-75'}, {'error': '1e-30'}, where)
        with self.assertRaises(ValueError):
            compare_records({'N': 'H^32'}, {'N': 'H^33'}, 'root.protocol')

    def test_existing_record_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'keep.json'; path.write_text('keep')
            with self.assertRaises(SystemExit):
                main(['--output', str(path)])
            self.assertEqual(path.read_text(), 'keep')


if __name__ == '__main__':
    unittest.main()
