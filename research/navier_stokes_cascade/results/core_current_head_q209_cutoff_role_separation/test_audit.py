import json
import tempfile
import unittest
from pathlib import Path
import audit
import q209


class Q209Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q209.statuses().values()))

    def test_source_role_separation(self):
        row = q209.source_role_audit()
        self.assertTrue(row["stage_specification_is_qualitative_smooth_plateau"])
        self.assertTrue(row["stage_field_formula_keeps_symbolic_chi"])
        self.assertFalse(row["stage_specification_names_squared_flat_formula"])
        self.assertTrue(row["older_diagnostic_uses_linear_reciprocal_exponent"])
        self.assertFalse(row["older_diagnostic_is_actual_supported_giant_field"])
        self.assertTrue(row["background_source_uses_squared_flat_step"])
        self.assertTrue(row["Q120_names_squared_flat_background_cutoff"])
        self.assertFalse(row["same_symbol_implies_same_function"])

    def test_distinct_cutoff_roles(self):
        row = q209.role_separation()
        self.assertIn("coefficient n", row["indices_are_distinct"])
        self.assertFalse(row["Q120_bit_budget_applies_to_Q205_stage_cutoff"])
        self.assertFalse(row["identifying_the_two_cutoffs_is_required_by_Q205"])
        self.assertFalse(row["historical_cutoff_changed_or_reselected"])

    def test_stage_cutoff_finiteness(self):
        row = q209.stage_cutoff_finiteness()
        self.assertEqual(row["order"], 36479350301)
        self.assertIn("S2", row["Euler_identity"])
        self.assertTrue(row["finite_sum_at_each_fixed_order"])
        self.assertTrue(row["compact_transition_and_continuity_give_finite_supremum"])
        self.assertFalse(row["numerical_K_order_evaluated"])
        self.assertFalse(row["Q120_formula_used"])

    def test_scale_absorption(self):
        row = q209.scale_absorption()
        self.assertTrue(row["chosen_after_current_constants"])
        self.assertTrue(row["endpoint_bound_holds"])
        self.assertTrue(row["decreasing_after_cutoff"])
        self.assertTrue(row["works_for_arbitrary_finite_stage_constants"])
        self.assertTrue(row["abstract_K_can_be_included_in_log_C"])
        self.assertFalse(row["stage_uniform_constant_used"])

    def test_corrected_frontier(self):
        row = q209.corrected_frontier()
        self.assertFalse(row["Q208_same_field_cutoff_identity_is_next_blocker"])
        self.assertEqual(row["R"], 36479350300)
        self.assertIn("T1/zeta", row["first_missing_majorant"])
        self.assertFalse(row["closed_A_R_gap_and_operator_bundle_recorded"])
        self.assertFalse(row["unconditional_residual_or_blowup"])

    def test_review_scope(self):
        row = q209.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["stage_cutoff_fixed_order_finiteness_proved"])
        self.assertFalse(row["scope"]["numerical_stage_cutoff_high_jet_bound_proved"])
        self.assertFalse(row["scope"]["Q120_budget_applied_to_stage_cutoff"])
        self.assertFalse(row["scope"]["complete_stage_602_primitive_bundle_proved"])

    def test_invalid_inputs(self):
        for call in (q209.statuses, q209.source_role_audit, q209.role_separation,
                     q209.stage_cutoff_finiteness, q209.scale_absorption,
                     q209.corrected_frontier, q209.review, q209.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q209.stage_cutoff_finiteness(-1)
        with self.assertRaises(ValueError):
            q209.stage_cutoff_finiteness(True)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 459)
        self.assertGreaterEqual(record["closure"]["edges"], 1677)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6168)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
