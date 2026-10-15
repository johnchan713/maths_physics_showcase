#!/usr/bin/env python3
"""Cutoff identities, exact summation budgets, failure controls and frozen scope."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from realization_bounds import (all_stage_order_review, comparison_budget, core_energy_review,
                                cutoff_loss, finite_field_loss, finite_residual_loss,
                                gain, increment_orders, scale_choice_review, tuple_loss)
from cutoff_checks import (cutoff_geometry_review, curl_review, divergent_series_review,
                           residual_review)
from audit import (FLAGS, HERE, compare_records, expected_protocol, generate,
                   inherited_provenance, main, validate)


class ExactBudgetTests(unittest.TestCase):
    def test_cycle_index_is_previous_B_not_next_B(self):
        first = increment_orders(1)
        self.assertEqual(first['wave_velocity'], F(7, 10)-F(1, 100000))
        self.assertEqual(first['mean_potential'], F(6, 5)-F(1, 50000))

    def test_all_stage_margins_have_equal_slopes(self):
        review = all_stage_order_review()
        self.assertEqual(review['row_slope'], review['target_slope'])
        self.assertGreater(min(review['all_stage_margins'].values()), 0)
        for j in (2, 101, 10**6):
            self.assertEqual({k: v-F(j, 10) for k, v in increment_orders(j).items()},
                             review['all_stage_margins'])

    def test_loss_reserves_one_full_curl_and_log_absorption(self):
        for m in (0, 2, 19):
            self.assertEqual(cutoff_loss(m), tuple_loss(m+1))
            self.assertEqual(finite_field_loss(m), cutoff_loss(m)+1)
            self.assertEqual(finite_residual_loss(m), F(251, 100)+2*m)

    def test_cutoff_tail_spends_half_the_increment_gain(self):
        self.assertEqual(gain(7)/2, F(7, 4000))

    def test_fixed_comparison_stage_covers_every_NS_term(self):
        for h in (F(1, 200), F(9, 1000), F(1, 10000)):
            for m, N in ((0, 0), (1, 3), (5, 11), (20, 2)):
                review = comparison_budget(m, N, h)
                self.assertTrue(all(review['checks'].values()))
                self.assertFalse(review['actual_stage_coefficients_or_scales_evaluated'])

    def test_conservative_example_is_not_actual_stage_evaluation(self):
        row = comparison_budget()
        self.assertEqual(row['comparison_stage'], 64079)
        self.assertEqual(row['old_field_times_tail_exponent'], 1)
        self.assertEqual(row['finite_residual_exponent'], F(29533, 1000))

    def test_low_stage_and_skipping_old_new_products_are_insufficient(self):
        L, H = cutoff_loss(2), finite_field_loss(2)
        self.assertLess(gain(2)/2-L-H, 0)
        # A tail barely decaying can still grow after multiplication by u_J.
        self.assertLess(F(1, 10)-H, 0)

    def test_invalid_stage_h_or_derivative_orders_are_rejected(self):
        for kwargs in ({'m': -1}, {'target': -1}, {'m': 1.5}, {'h': 0}, {'h': F(1, 100)}):
            with self.assertRaises(ValueError):
                comparison_budget(**kwargs)
        with self.assertRaises(ValueError):
            increment_orders(0)
        with self.assertRaises(ValueError):
            cutoff_loss(True)

    def test_all_interval_recipe_handles_extreme_synthetic_constants(self):
        review = scale_choice_review()
        self.assertTrue(all(review['checks'].values()))
        self.assertGreater(review['rows'][-1]['log_scale'], 10**14)
        self.assertFalse(review['actual_NS_scales_evaluated'])

    def test_energy_integrability_has_strict_uniform_h_margin(self):
        for h in (F(1, 200), F(999, 100000)):
            row = core_energy_review(h)
            self.assertTrue(all(row['checks'].values()))
            self.assertLess(row['axial_integrability_exponent'], 1)
            self.assertEqual(row['axial_slab_energy_power'], F(1, 2)-3*h)

    def test_raw_infinite_series_can_diverge_at_every_fixed_q(self):
        review = divergent_series_review()
        self.assertTrue(all(review['checks'].values()))
        self.assertFalse(FLAGS['raw_uncut_infinite_series_convergence_verified'])


class CutoffIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.geometry = cutoff_geometry_review()
        cls.curl = curl_review()
        cls.residual = residual_review()

    def test_actual_implicit_q_relation_holds_in_all_third_order_jets(self):
        self.assertTrue(self.geometry['actual_implicit_q'])
        self.assertTrue(self.geometry['checks']['implicit_q_equation_in_all_jet_coefficients'])

    def test_weighted_cutoff_jets_are_independent_of_scale(self):
        self.assertTrue(self.geometry['checks']['scale_independent_weighted_cutoff_jets'])

    def test_full_potential_product_rule_includes_two_velocity_derivatives(self):
        self.assertTrue(self.curl['checks']['full_curl_product_rule_through_two_velocity_derivatives'])

    def test_cut_wave_plus_mean_velocity_is_divergence_free(self):
        self.assertTrue(self.curl['checks']['complete_cut_velocity_is_divergence_free'])

    def test_independent_cartesian_and_cylindrical_evaluation(self):
        self.assertTrue(self.curl['checks']['cartesian_and_cylindrical_values_agree'])
        self.assertTrue(self.residual['checks']['independent_cartesian_residual_agrees'])

    def test_mean_radial_cutoff_commutator_has_correct_sign(self):
        self.assertTrue(self.curl['checks']['mean_cutoff_adds_negative_axial_derivative_radial_velocity'])
        self.assertTrue(self.curl['checks']['wrong_commutator_sign_is_detected'])

    def test_cutting_velocity_directly_fails(self):
        self.assertTrue(self.curl['checks']['naive_velocity_cutoff_breaks_divergence'])
        self.assertTrue(self.curl['checks']['curl_commutator_is_nonzero'])

    def test_both_cutoff_plateaus_and_zero_extensions_are_exact(self):
        self.assertTrue(self.curl['checks']['both_plateaus_are_exact'])
        self.assertTrue(self.curl['checks']['core_axis_and_exterior_support_preserved'])

    def test_complete_residual_increment_with_two_cutoff_scales(self):
        self.assertTrue(self.residual['two_distinct_cutoff_scales'])
        self.assertTrue(self.residual['checks']['exact_complete_residual_difference'])
        self.assertEqual(set(self.residual['retained_term_sizes']),
                         {'time', 'viscosity', 'pressure', 'old_new', 'new_old', 'new_square'})

    def test_all_interactions_are_nonzero(self):
        self.assertTrue(self.residual['checks']['every_linear_and_nonlinear_term_is_nonzero'])

    def test_cutoff_time_terms_are_not_omitted(self):
        self.assertTrue(self.residual['checks']['cutoff_time_terms_cannot_be_dropped'])

    def test_cutoff_viscosity_and_pressure_gradient_are_not_omitted(self):
        self.assertTrue(self.residual['checks']['cutoff_viscosity_commutators_cannot_be_dropped'])
        self.assertTrue(self.residual['checks']['pressure_cutoff_gradient_cannot_be_dropped'])


class EvidenceTests(unittest.TestCase):
    def test_protocol_and_all_inherited_source_bytes(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())
        inherited = inherited_provenance()
        self.assertGreater(len(inherited), 360)
        self.assertTrue(all(inherited.values()))

    def test_record_preserves_conditional_scope_and_open_global_gates(self):
        record = json.loads((HERE/'evidence.json').read_text())
        validate(record)
        self.assertEqual(len(record['gates']), 10)
        self.assertFalse(record['flags']['admissible_smooth_compact_force_verified'])
        self.assertFalse(record['flags']['blowup_verified'])
        self.assertTrue(record['flags']['conditional_local_kinetic_energy_bound_supplied'])

    def test_unconditional_or_forcing_promotion_is_rejected(self):
        for flag in ('blowup_verified', 'admissible_smooth_compact_force_verified',
                     'inherited_foundation_and_global_auxiliary_geometry_independently_verified'):
            record = json.loads((HERE/'evidence.json').read_text())
            record['flags'][flag] = True
            with self.assertRaises(ValueError):
                validate(record)

    def test_exact_quantifier_or_budget_tampering_is_rejected(self):
        record = json.loads((HERE/'evidence.json').read_text())
        bad = copy.deepcopy(record)
        bad['reviews']['fixed_stage_flat_residual_comparison']['rows'][0]['comparison_stage'] -= 1
        with self.assertRaises(ValueError):
            compare_records(record, bad)
        bad = copy.deepcopy(record)
        bad['conditional_proof']['comparison'] = 'let J grow with q'
        with self.assertRaises(ValueError):
            validate(bad)

    def test_frozen_replay_and_refusal_to_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'replay.json'
            main(['--verify-record', str(HERE/'evidence.json'), '--output', str(output)])
            self.assertTrue(output.exists())
            with self.assertRaises(SystemExit):
                main(['--output', str(output)])


if __name__ == '__main__':
    unittest.main()
