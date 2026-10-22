#!/usr/bin/env python3
"""Complete-primary norm routing and the initial-only Haar cancellation tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import covariance_majorant,fixture_budget,graph_majorant,review


class InitialScalarSeedBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_complete_real_covariance_norm_factor(self):
        self.assertEqual(covariance_majorant(0,Q(3),Q(4)),Q(25,2))
        self.assertEqual(covariance_majorant(3,Q(3),Q(4)),100)
        self.controls("primary_to_seed_budgets","complete_real_curl_square_reduction_margin")

    def test_graph_norm_includes_radial_and_connection_multipliers(self):
        self.assertEqual(graph_majorant(2,1,3,4),52)
        self.assertEqual(graph_majorant(2,2,3,4),68)

    def test_initial_pressure_and_axial_defects_use_complete_covariance(self):
        self.controls("initial_pressure_channels","complete_initial_pressure_scalar_identity")
        self.controls("initial_pressure_channels","complete_initial_axial_scalar_identity")
        self.controls("initial_pressure_channels","complete_initial_angular_scalar_nonzero")

    def test_compensator_changes_Haar_shape_but_not_zero_Haar_fast_input(self):
        self.controls("initial_pressure_channels","compensator_changes_full_Haar_pressure_shape")
        self.controls("initial_pressure_channels","compensator_does_not_change_zero_Haar_fast_pressure")
        self.controls("initial_pressure_channels","slow_stress_removed_only_by_Haar_projection")

    def test_compensator_and_pressure_compactness_are_retained(self):
        self.controls("initial_pressure_channels","both_moving_compensators_unit_mass")
        self.controls("initial_pressure_channels","both_full_Haar_pressures_compact")

    def test_whole_moving_pressure_moment_cannot_be_projected_out(self):
        self.controls("initial_pressure_channels","whole_moving_Haar_pressure_moment")
        self.controls("initial_pressure_channels","whole_moving_Haar_moment_derivative")
        self.controls("initial_pressure_channels","moving_compensator_term_cannot_be_deleted")
        self.controls("initial_pressure_channels","moving_cutoff_stays_in_zero_Haar_pressure_derivative")

    def test_initial_fast_zero_Haar_does_not_remove_its_new_squares(self):
        self.controls("fresh_initial_prefix","initial_fast_has_zero_Haar_all_components")
        self.controls("fresh_initial_prefix","actual_fast_quadratic_Haar_change_retained")

    def test_old_fast_slow_crosses_only_vanish_after_Haar(self):
        self.controls("fresh_initial_prefix","old_fast_slow_crosses_vanish_after_Haar_only")
        self.controls("fresh_initial_prefix","initial_sharper_Haar_velocities_are_slow_increment")

    def test_current_slow_remainders_are_independently_recomputed(self):
        self.controls("fresh_initial_prefix","fresh_slow_current_radial_remainder_exact")
        self.controls("fresh_initial_prefix","fresh_slow_angular_square_after_linear_row")
        self.controls("fresh_initial_prefix","fresh_slow_axial_square_after_linear_row")

    def test_initial_seed_defects_are_generally_nonzero(self):
        self.controls("fresh_initial_prefix","nonlinear_seed_defects_not_zero")
        self.controls("fresh_initial_prefix","slow_divergence_and_velocity_moments")

    def test_initial_simplification_fails_for_later_old_Haar_velocities(self):
        self.controls("fresh_initial_prefix","later_nonzero_old_Haar_breaks_initial_reduction")

    def test_positive_polynomial_chain_reaches_all_three_seed_scalars(self):
        self.controls("primary_to_seed_budgets","complete_primary_to_initial_scalar_polynomial_chain")
        self.controls("primary_to_seed_budgets","initial_Haar_and_full_velocity_majorants")

    def test_whole_axial_moment_uses_the_same_seed_polynomials(self):
        self.controls("primary_to_seed_budgets","whole_axial_moment_polynomial_chain")

    def test_C1_initial_budget_dominates_unreduced_dynamic_terms(self):
        self.controls("primary_to_seed_budgets","C1_initial_budget_dominates_unreduced_terms")
        self.controls("primary_to_seed_budgets","initial_C1_slow_square_ties_required_S2")

    def test_initial_fast_and_Haar_rates_are_distinct(self):
        self.controls("primary_to_seed_budgets","initial_fast_squares_still_supply_S1")
        self.controls("primary_to_seed_budgets","Haar_rates_are_sharper_than_full_fast_rates")
        self.controls("primary_to_seed_budgets","initial_radial_time_and_base_products_are_higher_order")

    def test_primary_coefficient_input_orders_are_traced_through_actual_calls(self):
        self.controls("primary_to_seed_budgets","scalar_primary_input_order_is_n_plus_10")
        self.controls("primary_to_seed_budgets","axial_moment_primary_input_order_is_n_plus_11")
        self.controls("primary_to_seed_budgets","input_order_ledger")

    def test_aggregate_label_bound_is_not_inferred_from_individual_bounds(self):
        self.controls("primary_to_seed_budgets","aggregate_label_norm_not_inferred_from_individual_bounds")

    def test_invalid_orders_geometry_norms_and_scale_are_rejected(self):
        with self.assertRaises(ValueError):covariance_majorant(-1,1,1)
        with self.assertRaises(ValueError):graph_majorant(1,-1,1,1)
        with self.assertRaises(ValueError):graph_majorant(1,1,-1,1)
        with self.assertRaises(TypeError):graph_majorant(1,1,.1,1)
        with self.assertRaises(ValueError):fixture_budget(0).output(0)
        with self.assertRaises(ValueError):fixture_budget().output(True)

    def test_both_pressure_memories_and_complete_covariance_remain(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))
        self.controls("complete_covariance","all_six_covariances_nonzero")
        self.controls("complete_covariance","omitted_axial_velocity_detected")

    def test_frozen_historical_sources_and_new_record_reproduce(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_initial_prefix_is_not_promoted_to_full_canonical_seed(self):
        self.assertTrue(self.record["flags"]["known_qualitative_initial_S2_order_not_presented_as_new"])
        for key in ("primary_aggregate_root_and_edge_hypotheses_independently_verified","full_signed_sourced_canonical_finite_cycle_seed_verified","complete_pointwise_mean_or_full_residual_verified","uniform_infinite_stage_or_all_order_constant_proved","actual_profile_tilt_width_or_accepted_tail_changed","unconditional_pulse_or_blowup_proved"):
            self.assertFalse(self.record["flags"][key])
        changed=deepcopy(self.record);changed["flags"]["full_signed_sourced_canonical_finite_cycle_seed_verified"]=True
        with self.assertRaises(ValueError):audit.validate(changed)

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"evidence.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
