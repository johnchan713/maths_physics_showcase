#!/usr/bin/env python3
"""Compact mean equations, exact increments, inverse domains and scope."""
import copy
from fractions import Fraction
import json
import tempfile
from pathlib import Path
import unittest
import mpmath as mp
import numpy as np
from operators import (divisor, mean_update_budget, moment_blocks, radial_mode,
                       solve_moment_targets, temporal_inverse)
from reviews import (defect_increment_review, integrated_mean_review, moment_map_review,
                     physical_mean_review, torus_and_radial_review)
from audit import (FLAGS, HERE, compare_records, expected_protocol,
                   inherited_provenance, main, validate)


class CompleteMeanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.review = physical_mean_review()
        cls.defects = defect_increment_review()

    def test_full_curl_covariance_mean_matches_physical_momentum(self):
        self.assertTrue(self.review['checks']['full_six_entry_covariance_mean_identity'])

    def test_normalized_time_axial_and_viscosity_factors(self):
        self.assertTrue(self.review['checks']['normalized_time_axial_and_viscous_factors'])

    def test_every_remaining_axial_covariance_entry_is_used(self):
        for key in ('axial_cross_covariance', 'axial_normal_covariance', 'radial_axial_covariance'):
            self.assertGreater(self.review['omission_gaps'][key], 1e-5)

    def test_mean_velocity_is_an_exact_compact_curl(self):
        self.assertTrue(self.review['checks']['mean_curl_is_divergence_free'])

    def test_radial_memory_viscosity_and_new_square_are_not_dropped(self):
        for key in ('radial_background_memory', 'viscosity', 'quadratic_increment'):
            self.assertGreater(self.review['omission_gaps'][key], 1e-5)

    def test_exact_radial_increment_identity(self):
        self.assertTrue(self.review['checks']['exact_increment_retains_radial_memory_and_quadratic_terms'])

    def test_all_three_defects_are_recomputed_from_same_field(self):
        self.assertTrue(self.defects['checks']['all_three_exact_defect_increment_identities'])

    def test_linear_target_solve_does_not_remove_nonlinear_defects(self):
        self.assertTrue(self.defects['checks']['nonlinear_defects_do_not_vanish_after_a_linear_solve'])


class CompactInverseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.review = torus_and_radial_review()

    def test_actual_torus_eigendirections_and_algebraic_divisor_bound(self):
        self.assertTrue(self.review['checks']['actual_auxiliary_eigendirections'])
        self.assertTrue(self.review['checks']['algebraic_divisor_bound_crosscheck'])

    def test_temporal_inverse_cancels_nonzero_modes_and_keeps_reality(self):
        self.assertTrue(self.review['checks']['zero_mean_temporal_inverse_and_reality'])

    def test_temporal_inverse_agrees_on_actual_integer_coverings(self):
        matrix = mp.matrix([[3, 1], [1, 5]])
        tg = 4+mp.sqrt(2)
        for index in (0, 2, 5, 8):
            cover = matrix**index
            for mode in ((1, -2), (-3, 1), (2, 1)):
                lifted = tuple(int(x) for x in cover.T*mp.matrix(mode))
                self.assertLess(abs(divisor(lifted, True)/tg**index-divisor(mode, True)), mp.mpf('1e-70'))
                physical = temporal_inverse({lifted: 1}, 1)[lifted]
                chart = tg**(-index)*temporal_inverse({mode: 1}, 1)[mode]
                self.assertLess(abs(physical-chart), mp.mpf('1e-70'))

    def test_auxiliary_zero_mode_must_be_removed_first(self):
        with self.assertRaises(ValueError):
            temporal_inverse({(0, 0): 1, (1, 0): 2}, 1)

    def test_invalid_temporal_scale_is_not_accepted(self):
        for ci in (0, -1, 'inf', 'nan'):
            with self.assertRaises(ValueError):
                temporal_inverse({(1, 0): 1}, ci)

    def test_all_cylindrical_weighted_radial_graph_identities(self):
        self.assertTrue(self.review['checks']['all_three_weighted_radial_graph_identities'])

    def test_phase_path_omission_fails(self):
        self.assertTrue(self.review['checks']['radial_path_cannot_be_omitted'])

    def test_radial_compact_support_values(self):
        self.assertTrue(self.review['checks']['primitive_values_keep_shell_support'])

    def test_repeated_integration_by_parts_keeps_the_exact_full_integral(self):
        self.assertTrue(self.review['checks']['integration_by_parts_identity'])

    def test_nonzero_weighted_zero_mode_is_a_real_obstruction(self):
        self.assertTrue(self.review['checks']['zero_radial_mean_condition_is_necessary'])

    def test_invalid_radial_weight_is_rejected(self):
        with self.assertRaises(ValueError):
            radial_mode(1.4, (1, 2), 8, 3, lambda r: .5)


class PressureAndMomentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.integrated = integrated_mean_review()
        cls.blocks = moment_map_review()

    def test_pressure_source_has_zero_integral_and_compact_pressure_moment(self):
        for key in ('compact_pressure_source_has_zero_integral', 'pressure_moment_keeps_radial_source'):
            self.assertTrue(self.integrated['checks'][key])

    def test_complete_integrated_axial_equation(self):
        self.assertTrue(self.integrated['checks']['complete_integrated_mean_identities'])

    def test_moving_pressure_bump_and_pressure_axial_moment_are_needed(self):
        self.assertTrue(self.integrated['checks']['moving_pressure_bump_derivative_required'])
        self.assertTrue(self.integrated['checks']['axial_pressure_moment_required'])

    def test_both_velocity_moments_vanish(self):
        self.assertTrue(self.integrated['checks']['two_velocity_moments_are_preserved'])

    def test_without_zero_moment_time_and_axial_viscosity_terms_survive(self):
        self.assertTrue(self.integrated['checks']['arbitrary_nonzero_velocity_moment_cannot_be_dropped'])
        self.assertTrue(self.integrated['checks']['off_constraint_time_and_axial_viscosity_identity'])

    def test_refined_stencil_and_edge_split_quadrature_reproduce(self):
        refined = integrated_mean_review(nodes=128, step=.001)
        self.assertTrue(all(refined['checks'].values()))
        for key in ('integrated_mean_residual', 'expected_integrated_flux_derivatives',
                    'radial_defect', 'pressure_moment', 'c_rho_z'):
            np.testing.assert_allclose(self.integrated[key], refined[key], rtol=1e-7, atol=2e-9)

    def test_all_five_signed_target_rows_and_scaled_bump_determinants(self):
        self.assertTrue(self.blocks['checks']['all_five_signed_target_rows'])
        self.assertTrue(self.blocks['checks']['exact_scaled_bump_determinants'])

    def test_large_signed_targets_do_not_require_a_nonlinear_root(self):
        at, az = moment_blocks('.2', '1.3')
        u, s = solve_moment_targets(mp.mpf('1e30'), -mp.mpf('2e30'), mp.mpf('3e30'))
        self.assertLess(mp.norm(at*u-mp.matrix([0, -mp.mpf('1e30'), -mp.mpf('3e30')]))/mp.mpf('1e30'), mp.mpf('1e-60'))
        self.assertLess(mp.norm(az*s-mp.matrix([0, mp.mpf('2e30')]))/mp.mpf('1e30'), mp.mpf('1e-60'))

    def test_lambda_zero_or_patch_amplitude_zero_is_not_invertible(self):
        for lam, amplitude in ((0, 1), (-1, 1), ('.2', 0)):
            with self.assertRaises(ValueError):
                moment_blocks(lam, amplitude)
        self.assertTrue(self.blocks['checks']['fixed_small_lambda_has_costly_inverse'])

    def test_mean_gain_retains_every_radial_source_term(self):
        budget = mean_update_budget()
        self.assertEqual(len(budget['radial_source_terms']), 11)
        self.assertEqual(budget['moment_defect_gain'], Fraction(44999, 50000))
        self.assertGreaterEqual(budget['minimum_margin'], 0)
        self.assertFalse(budget['full_residual_cycle_gain_verified'])

    def test_mean_gain_rejects_missing_current_size_hypothesis(self):
        with self.assertRaises(ValueError):
            mean_update_budget(alpha=Fraction(1, 2))


class EvidenceTests(unittest.TestCase):
    def test_frozen_protocol_and_inherited_source_bytes(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())
        provenance = inherited_provenance()
        self.assertGreater(len(provenance), 100)
        self.assertTrue(all(provenance.values()))

    def test_frozen_record_does_not_assert_global_construction(self):
        record = json.loads((HERE/'evidence.json').read_text())
        validate(record)
        self.assertFalse(FLAGS['full_residual_cycle_gain_verified'])
        self.assertFalse(FLAGS['blowup_verified'])

    def test_claim_promotion_or_missing_gate_is_rejected(self):
        record = json.loads((HERE/'evidence.json').read_text())
        for change in ('blowup', 'cycle', 'gate'):
            wrong = copy.deepcopy(record)
            if change == 'gate':
                wrong['gates'].pop('same_field_nonlinear_defect_changes')
            else:
                wrong['flags']['blowup_verified' if change == 'blowup' else 'full_residual_cycle_gain_verified'] = True
            with self.assertRaises(ValueError):
                validate(wrong)

    def test_frozen_numerical_or_source_mutation_is_rejected(self):
        record = json.loads((HERE/'evidence.json').read_text())
        for field in ('source_hashes', 'reviews'):
            wrong = copy.deepcopy(record)
            if field == 'source_hashes':
                key = next(iter(wrong[field]))
                wrong[field][key] = '0'*64
            else:
                wrong[field]['pressure_and_integrated_constraints']['radial_defect'] += 1
            with self.assertRaises(ValueError):
                compare_records(record, wrong)

    def test_existing_output_is_protected(self):
        with tempfile.TemporaryDirectory(dir=HERE) as folder:
            path = Path(folder)/'already-exists.json'
            path.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output', str(path)])
            self.assertEqual(path.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
