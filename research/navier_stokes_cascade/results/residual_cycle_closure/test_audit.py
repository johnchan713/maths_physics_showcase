#!/usr/bin/env python3
"""Nonlinear identities, all-stage order inequalities and conditional scope."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from cycle_bounds import (Affine, KAPPA, MIN_B, STEP, cycle_certificate,
                          initialization, line, physical_loss, physical_residual_power)
from interaction_checks import (covariance_increment_review, full_increment_review,
                                longitudinal_and_transport_review)
from support_checks import sharp_loss_controls, support_and_flat_review
from audit import (FLAGS, HERE, compare_records, expected_protocol,
                   inherited_provenance, main, validate)


class ExactCycleTests(unittest.TestCase):
    def test_initial_state_meets_wave_mean_defect_and_total_bounds(self):
        result = initialization()
        self.assertTrue(all(v > 0 for v in result['margins'].values()))
        self.assertEqual(result['orders']['updated_tangential_mean'], 2-3*KAPPA)
        self.assertEqual(result['orders']['updated_compatibility_defects'], F(19, 10)-3*KAPPA)

    def test_cycle_inequalities_hold_on_the_entire_stage_half_line(self):
        result = cycle_certificate()
        self.assertTrue(all(v >= 0 for v in result['all_half_line_margins'].values()))
        self.assertEqual(result['minimum_margin'], F(499, 50000))
        self.assertEqual(result['requested_gain'], F(1, 10))

    def test_large_B_has_no_hidden_negative_slope(self):
        result = cycle_certificate()
        for name, rows in result['term_groups'].items():
            target = (line(1, F(67, 100)) if name == 'averaged_covariance_remainders'
                      else result['mean_and_defect_target'] if name == 'final_mean_rows'
                      else result['wave_target'])
            for value in rows.values():
                self.assertGreaterEqual(value.slope-target.slope, 0)
                self.assertGreaterEqual(value.margin_above(target), 0)

    def test_finite_stage_success_does_not_certify_a_bad_half_line(self):
        bad = line(0, 100)
        target = line(1, 0)
        self.assertGreater(bad.at(MIN_B), target.at(MIN_B))
        with self.assertRaises(ValueError):
            bad.margin_above(target)

    def test_signed_square_and_all_interaction_rows_are_present(self):
        groups = cycle_certificate()['term_groups']
        self.assertEqual(set(groups['step1_wave']), {'linear', 'old_wave_cross', 'new_wave_square', 'total_mean_cross'})
        self.assertEqual(set(groups['step2_wave']), {'linear', 'total_mean_cross', 'old_wave_cross', 'signed_square'})
        self.assertIn('signed_square', groups['averaged_covariance_remainders'])
        self.assertIn('axial_pressure_and_flux', groups['averaged_covariance_remainders'])

    def test_unproved_larger_gain_or_large_graph_loss_is_rejected(self):
        for kwargs in ({'requested_gain': F(1, 5)}, {'kappa': F(1, 100)}, {'kappa': 0}):
            with self.assertRaises(ValueError):
                cycle_certificate(**kwargs)

    def test_total_mean_and_fixed_primary_remainder_bounds_are_preserved(self):
        margins = cycle_certificate()['all_half_line_margins']
        self.assertGreater(margins['preserve_total_wave_remainder'], 0)
        self.assertGreater(margins['preserve_total_mean'], 0)
        self.assertGreater(margins['five_map_input_at_least_0.9'], 0)

    def test_physical_order_gain_survives_a_fixed_derivative_loss(self):
        for h in (F(1, 200), F(1, 500)):
            for m in (0, 1, 3, 8):
                self.assertEqual(physical_loss(h, m), F(5, 2)+2*h+2*m)
                for j in (0, 3, 1000):
                    self.assertEqual(physical_residual_power(j+1, m, h)-physical_residual_power(j, m, h), h*STEP)

    def test_initial_physical_bound_is_not_called_absolute_decay(self):
        self.assertEqual(physical_residual_power(0), -F(5013, 2000))
        self.assertLess(physical_residual_power(100), 0)
        self.assertFalse(FLAGS['low_stage_physical_residual_proved_to_decay'])
        self.assertFalse(FLAGS['fixed_scale_residual_norm_contraction_verified'])


class NonlinearIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.transport = longitudinal_and_transport_review()
        cls.full = full_increment_review()
        cls.covariance = covariance_increment_review()

    def test_longitudinal_identity_uses_complete_curl_velocity(self):
        self.assertTrue(self.transport['checks']['full_curl_longitudinal_identity'])
        self.assertTrue(self.transport['checks']['exact_complex_curl_incompressibility'])

    def test_transport_matches_independent_cartesian_curls(self):
        self.assertTrue(self.transport['checks']['direct_cartesian_transport_agrees'])
        self.assertTrue(self.transport['checks']['all_harmonic_signs_and_zero_output_checked'])

    def test_principal_only_velocity_and_missing_frame_terms_fail(self):
        self.assertTrue(self.transport['checks']['dropping_curl_breaks_divergence'])
        self.assertTrue(self.transport['checks']['cylindrical_connection_is_present'])

    def test_same_label_fast_loss_is_canceled_and_uncontrolled_longitudinal_part_fails(self):
        self.assertTrue(self.transport['checks']['same_label_transport_has_no_half_power_loss'])
        self.assertTrue(self.transport['checks']['without_longitudinal_control_fast_loss_returns'])

    def test_full_same_field_residual_increment(self):
        self.assertTrue(self.full['checks']['full_actual_velocity_increment_identity'])
        self.assertTrue(self.full['checks']['nonlinear_increment_has_all_three_products'])

    def test_pressure_viscosity_old_new_and_new_square_are_nonzero(self):
        for key in ('new_velocity_square_is_not_discarded', 'pressure_is_not_discarded',
                    'viscosity_is_not_discarded', 'old_new_cross_terms_are_not_discarded'):
            self.assertTrue(self.full['checks'][key])

    def test_signed_primary_inverse_target(self):
        self.assertTrue(self.covariance['checks']['signed_primary_linear_inverse_matches_target'])

    def test_signed_covariance_keeps_every_actual_velocity_remainder(self):
        for key in ('complete_covariance_change_has_all_remainders', 'old_wave_remainder_is_retained',
                    'new_curl_remainder_is_retained', 'signed_increment_square_is_retained'):
            self.assertTrue(self.covariance['checks'][key])


class SourceAndHypothesisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.support = support_and_flat_review()
        cls.sharp = sharp_loss_controls()

    def test_temporal_inverse_can_fill_the_auxiliary_torus(self):
        self.assertTrue(self.support['checks']['temporal_inverse_can_fill_auxiliary_torus'])
        self.assertTrue(self.support['checks']['temporal_inverse_fourier_identity'])

    def test_mean_wave_product_still_has_wave_support(self):
        self.assertTrue(self.support['checks']['mean_wave_product_keeps_wave_support'])

    def test_separated_closed_supports_have_zero_products(self):
        self.assertTrue(self.support['checks']['separated_support_products_vanish'])

    def test_mean_zero_flat_coefficient_must_stay_in_later_products(self):
        self.assertTrue(self.support['checks']['zero_mean_flat_velocity_is_retained_in_later_products'])

    def test_finite_angular_harmonics_do_not_depend_on_band(self):
        self.assertTrue(self.support['checks']['angular_harmonics_remain_finite_and_band_independent'])

    def test_overlapping_different_phases_restore_the_fast_loss(self):
        self.assertTrue(self.sharp['checks']['closed_support_separation_is_needed'])
        self.assertTrue(self.sharp['checks']['same_label_incompressibility_cannot_be_used_for_different_phases'])

    def test_mean_advection_half_power_loss_cannot_be_ignored(self):
        self.assertTrue(self.sharp['checks']['mean_wave_half_power_loss_is_sharp'])
        self.assertTrue(self.sharp['checks']['dropping_mean_advection_of_phase_is_invalid'])


class EvidenceTests(unittest.TestCase):
    def test_protocol_and_all_inherited_source_hashes(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())
        inherited = inherited_provenance()
        self.assertGreater(len(inherited), 350)
        self.assertTrue(all(inherited.values()))

    def test_evidence_states_the_remaining_realization_and_foundation_gaps(self):
        record = json.loads((HERE/'evidence.json').read_text())
        validate(record)
        self.assertFalse(record['flags']['complete_infinite_wave_and_mean_potential_sum_verified'])
        self.assertFalse(record['flags']['inherited_foundation_and_global_auxiliary_geometry_independently_verified'])
        self.assertFalse(record['flags']['blowup_verified'])

    def test_unconditional_claim_or_missing_gate_is_rejected(self):
        record = json.loads((HERE/'evidence.json').read_text())
        for key in ('unconditional_complete_residual_construction_verified', 'blowup_verified'):
            wrong = copy.deepcopy(record)
            wrong['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(wrong)
        wrong = copy.deepcopy(record)
        wrong['gates'].pop('source_support_and_flat_memory')
        with self.assertRaises(ValueError):
            validate(wrong)

    def test_frozen_exact_orders_or_numerical_observations_cannot_be_changed(self):
        record = json.loads((HERE/'evidence.json').read_text())
        for part in ('all_stage_cycle_orders', 'longitudinal_and_transport'):
            wrong = copy.deepcopy(record)
            if part == 'all_stage_cycle_orders':
                wrong['reviews'][part]['minimum_margin'] = '1/100'
            else:
                wrong['reviews'][part]['principal_only_divergence_gap'] += 1
            with self.assertRaises(ValueError):
                compare_records(record, wrong)

    def test_preexisting_output_is_preserved(self):
        with tempfile.TemporaryDirectory(dir=HERE) as folder:
            path = Path(folder)/'preserve.json'
            path.write_text('keep')
            with self.assertRaises(SystemExit):
                main(['--output', str(path)])
            self.assertEqual(path.read_text(), 'keep')


if __name__ == '__main__':
    unittest.main()
