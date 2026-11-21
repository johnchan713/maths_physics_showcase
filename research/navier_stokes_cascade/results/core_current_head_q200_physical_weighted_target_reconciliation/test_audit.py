import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q200


class Q200Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q200.statuses().values()))

    def test_corrected_stage_rates_match(self):
        row = q200.rates()
        self.assertEqual(row["Q197_current_target_input_C"], Q(149999, 100000))
        self.assertEqual(row["Q197_post_source_regular_D"], Q(74999, 50000))
        self.assertEqual(row["Q197_post_source_regular_D"], row["Q160_first_source_B0"])
        self.assertEqual(row["Q198_signed_output_C0"], Q(99999, 50000))
        self.assertEqual(row["Q198_signed_output_C0"], row["Q160_target_seed_C0"])

    def test_same_physical_stage_without_family_promotion(self):
        row = q200.stage_identification()
        self.assertTrue(row["Q197_D_equals_Q160_B0"])
        self.assertTrue(row["Q198_C0_equals_Q160_C0"])
        self.assertEqual(row["initial_target_is_stronger_than_Q197_C"], Q(24999, 50000))
        self.assertTrue(row["same_inherited_physical_family_only"])
        self.assertFalse(row["physical_and_canonical_families_identified"])

    def test_Q199_obstruction_is_retained_with_correct_scope(self):
        route = q200.actual_weight_route()
        self.assertTrue(route["Q199_scaled_bump_obstruction_retained"])
        self.assertFalse(route["standalone_unweighted_to_weighted_inference_valid"])
        self.assertFalse(route["raw_Q199_inequality_proved_as_generic_seminorm_implication"])
        self.assertTrue(route["actual_same_family_target_constructed_by_state_difference"])

    def test_direct_route_retains_weight_and_pressure_memory(self):
        route = q200.actual_weight_route()
        for key in (
            "direct_physical_route_present",
            "cancellation_before_estimation",
            "all_six_covariances_in_pressure_reconstruction",
            "axial_velocity_retained",
            "compact_pressure_retained",
            "whole_moving_pressure_product_retained",
            "moving_projection_commutator_retained",
            "reserved_bumps_use_positive_interior_weight_floor",
        ):
            self.assertTrue(route[key])

    def test_generic_signed_errors_tie_at_C0(self):
        row = q200.generic_step()
        self.assertTrue(row["all_errors_equal_output"])
        self.assertEqual(set(row["errors"].values()), {Q(99999, 50000)})
        self.assertEqual(row["gain"], Q(49999, 100000))

    def test_future_cycle_is_simultaneous_with_two_kappa_slack(self):
        for j in range(8):
            row = q200.simultaneous_step(j)
            self.assertTrue(row["source_and_target_advanced_simultaneously"])
            self.assertFalse(row["sequential_Q160_then_source_reading_valid"])
            self.assertEqual(row["target_slack"], Q(1, 50000))
            self.assertTrue(row["closes"])

    def test_Q199_global_conclusion_corrected(self):
        row = q200.correction()
        self.assertTrue(row["Q199_local_provider_cut_was_accurately_reported"])
        self.assertTrue(row["Q199_global_no_route_conclusion_was_too_strong"])
        self.assertIn("core_current_head_q160_physical_target_induction", row["omitted_physical_providers"])
        self.assertIn("core_current_head_q180_target_source_simultaneous_induction", row["omitted_physical_providers"])

    def test_physical_limit_and_frontier(self):
        self.assertEqual(q200.physical_residual_power(), Q(-15025001, 10000000))
        self.assertLess(q200.physical_residual_power(), 0)
        boundary = q200.frontier()
        self.assertTrue(boundary["same_family_weighted_target_route_available"])
        self.assertFalse(boundary["post_signed_means_recomputed_here"])
        self.assertFalse(boundary["complete_residual_or_blowup_proved_here"])

    def test_review(self):
        row = q200.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["analytical_relative_to_named_continuum_providers"])
        self.assertFalse(row["scope"]["numerical_or_symbolic_tests_are_continuum_proof"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q200.statuses, q200.stage_identification, q200.actual_weight_route, q200.correction, q200.frontier, q200.review, q200.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError): q200.rates(0)
        with self.assertRaises(ValueError): q200.simultaneous_step(-1)
        with self.assertRaises(ValueError): q200.physical_residual_power(Q(1, 50))

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 448)
        self.assertGreaterEqual(record["closure"]["edges"], 1602)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6048)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
