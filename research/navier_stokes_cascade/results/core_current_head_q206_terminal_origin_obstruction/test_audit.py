import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q206


class Q206Tests(unittest.TestCase):
    def test_provider_statuses_and_same_object(self):
        self.assertTrue(all(q206.statuses().values()))
        self.assertTrue(all(q206.same_object_links().values()))

    def test_protected_geometry(self):
        r = q206.geometry()
        self.assertEqual(r["h"], Q(1, 200))
        self.assertEqual(r["A"], Q(101, 200))
        self.assertEqual(r["two_A"], Q(101, 100))
        self.assertEqual(r["support_margin"], Q(3, 101))
        self.assertEqual(r["inverse_2h"], 100)
        self.assertTrue(r["annular_generated_support_starts_at_Xa"])
        self.assertTrue(r["localization_plateau_preserves_both_paths"])

    def test_velocity_value_jet_obstruction(self):
        r = q206.velocity_obstruction()
        self.assertEqual(r["threshold_power"], 100)
        self.assertEqual(r["lower_bound_exponent"], Q(-101, 200))
        self.assertTrue(r["path_tends_to_terminal_origin"])
        self.assertFalse(r["continuous_terminal_velocity_possible"])
        self.assertFalse(r["value_jet_exists"])

    def test_pressure_two_path_obstruction(self):
        r = q206.pressure_obstruction()
        self.assertEqual(r["X2"], Q(3, 202))
        self.assertEqual(r["log_ratio_lower"], Q(2, 5))
        self.assertEqual(r["d0_rational_lower"], Q(9, 500))
        self.assertEqual(r["lower_bound_exponent"], Q(-101, 100))
        self.assertTrue(r["time_dependent_spatial_gauge_cancels"])
        self.assertFalse(r["continuous_terminal_pressure_in_any_time_gauge_possible"])
        self.assertFalse(r["value_jet_exists"])

    def test_endpoint_force_separation(self):
        r = q206.endpoint_separation()
        self.assertTrue(r["punctured_one_sided_velocity_pressure_jets"])
        self.assertFalse(r["terminal_origin_in_punctured_atlas"])
        self.assertTrue(r["smooth_localized_force_trace"])
        self.assertFalse(r["smooth_force_implies_smooth_field"])
        self.assertTrue(r["complete_residual_corner_flatness"])
        self.assertFalse(r["corner_flatness_implies_field_value_bound"])
        self.assertFalse(r["future_velocity_or_pressure_constructed"])

    def test_review_scope(self):
        r = q206.review()
        self.assertTrue(all(r["checks"].values()))
        self.assertTrue(r["scope"]["regular_terminal_velocity_obstructed"])
        self.assertTrue(r["scope"]["regular_terminal_pressure_obstructed"])
        self.assertTrue(r["scope"]["punctured_terminal_jets_preserved"])
        self.assertTrue(r["scope"]["smooth_force_trace_preserved"])
        self.assertFalse(r["scope"]["unconditional_complete_residual_force_exclusion_or_blowup"])

    def test_invalid_identity_and_inputs(self):
        for call in (q206.statuses, q206.same_object_links, q206.geometry,
                     q206.velocity_obstruction, q206.pressure_obstruction,
                     q206.endpoint_separation, q206.review, q206.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q206.pressure_obstruction(e0=0)
        with self.assertRaises(ValueError):
            q206.pressure_obstruction(X2=Q(5, 101))

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 456)
        self.assertGreaterEqual(record["closure"]["edges"], 1655)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6133)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
