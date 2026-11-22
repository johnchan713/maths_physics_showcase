import json
import tempfile
import unittest
from pathlib import Path
import audit
import q208


class Q208Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q208.statuses().values()))

    def test_frozen_same_object_paths(self):
        rows = q208.same_object_links()
        self.assertTrue(all(value > 0 for value in rows.values()))
        self.assertEqual(rows["Q207_to_Q205"], 2)
        self.assertGreaterEqual(rows["Q207_to_Q120"], 1)

    def test_exact_existing_cutoff_budget(self):
        row = q208.cutoff_provider()
        self.assertEqual(row["identity"], "implemented-squared-flat-Borel")
        self.assertIn("exp(-1/x^2)", row["flat_seed"])
        self.assertEqual(row["ceil_log2_order"], 36)
        self.assertEqual(row["binary_budget"], 4414001386307)
        self.assertTrue(row["direct_q120_budget_match"])
        self.assertTrue(row["stage_602_order_match"])
        self.assertFalse(row["factorial_or_derivative_array_materialized"])
        self.assertFalse(row["new_cutoff_substituted"])
        self.assertFalse(row["q120_actual_field_identity_certified"])

    def test_correction_and_schedule(self):
        row = q208.correction()
        self.assertEqual(row["discarded_draft_seed"], "exp(-1/x)")
        self.assertTrue(row["Q120_already_proved_all_order_cutoff_bound"])
        self.assertTrue(row["Q204_rule_selects_background_scale_after_finite_constants"])
        self.assertTrue(row["Q205_rule_selects_stage_scale_after_finite_constants"])
        self.assertTrue(row["conditional_budget_available"])
        self.assertFalse(row["Q207_to_Q120_same_field_identity_explicitly_certified"])
        self.assertFalse(row["cutoff_shape_removed_from_frontier_unconditionally"])
        self.assertFalse(row["historical_evidence_modified"])

    def test_propagation_and_frontier(self):
        propagation = q208.propagation_ledger()
        frontier = q208.remaining_frontier()
        self.assertTrue(propagation["cutoff_shape_bound_available_after_identity_binding"])
        self.assertFalse(propagation["same_field_cutoff_identity_bound"])
        self.assertFalse(propagation["selected_profile_gap_and_operator_bundle_closed_directly"])
        self.assertEqual(frontier["R"], 36479350300)
        self.assertIn("T1/zeta", frontier["first_missing_majorant"])
        self.assertIn("compact/slow inverse", frontier["additional_required_bundle"])
        self.assertIn("bind Q207", frontier["same_cutoff_identity_obligation"])
        self.assertFalse(frontier["closed_majorant_A_R_from_selected_defining_maps_recorded"])
        self.assertFalse(frontier["unconditional_residual_or_blowup"])

    def test_review_scope(self):
        row = q208.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["existing_scale_cutoff_bound_available_conditionally"])
        self.assertEqual(row["scope"]["stage_602_cutoff_binary_budget"], 4414001386307)
        self.assertFalse(row["scope"]["same_field_cutoff_identity_certified"])
        self.assertFalse(row["scope"]["complete_stage_602_primitive_bundle_proved"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q208.statuses, q208.same_object_links, q208.cutoff_provider,
                     q208.correction, q208.propagation_ledger,
                     q208.remaining_frontier, q208.review, q208.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q208.cutoff_provider(2)
        with self.assertRaises(ValueError):
            q208.cutoff_provider(True)
        with self.assertRaises(ValueError):
            q208.frozen_path(q208.PROVIDERS[0], "missing-provider")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 458)
        self.assertGreaterEqual(record["closure"]["edges"], 1670)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6157)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
