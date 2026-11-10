import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q157


class Q157Tests(unittest.TestCase):
    def test_fixed_parameters(self):
        row = q157.fixed_parameters()
        self.assertEqual(row["h"], Q(1, 200))
        self.assertEqual(row["A"], Q(101, 200))
        self.assertEqual(row["X_a"], Q(4, 101))
        self.assertEqual(row["X_in"], Q(1, 101))
        self.assertLess(row["X_in"], row["X_a"])

    def test_exact_path(self):
        row = q157.path_coordinates(Q(1, 1000))
        self.assertEqual(row["q"], Q(1, 1000))
        self.assertEqual(row["eta"], 0)
        self.assertEqual(row["X"], Q(1, 101))
        self.assertEqual(row["r_squared"], Q(2, 101000))

    def test_core_threshold(self):
        row = q157.combined_threshold(
            Q(1, 10), 10, 1, Q(1, 2), Q(1, 4)
        )
        self.assertGreater(row["q_core"], 0)
        self.assertGreater(row["plateau"], 0)
        self.assertEqual(row["tau_star"], min(row["q_core"], row["plateau"]))

    def test_same_physical_binding(self):
        row = q157.same_physical_binding()
        self.assertTrue(row["path_inside_hole"])
        self.assertFalse(row["canonical_physical_primary_identification_used"])
        self.assertFalse(row["pressure_gauge_changed"])
        self.assertIn("Q153", row["diagonal"])

    def test_lower_bound_is_quantitative(self):
        row = q157.lower_bound_statement()
        self.assertEqual(row["exponent"], Q(101, 200))
        self.assertIn("tau_star", row["bound"])
        self.assertIn("(e0/(2 C_bg))^100", row["tau_star"])
        self.assertTrue(row["conditional"])

    def test_conditional_exclusion_scope(self):
        row = q157.exclusion_statement()
        self.assertTrue(row["global_smooth_finite_energy_competitor_excluded_conditionally"])
        self.assertFalse(row["unconditional_blowup_claimed"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])

    def test_provider_contracts(self):
        row = q157.provider_contracts()
        self.assertTrue(all(row.values()))

    def test_record_scope(self):
        row = q157.record()
        self.assertTrue(row["same_physical_core_and_conditional_exclusion_certified"])
        self.assertFalse(row["unconditional_admissible_force_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            q157.fixed_parameters("wrong")
        with self.assertRaises(ValueError):
            q157.path_coordinates(0)
        with self.assertRaises(ValueError):
            q157.core_threshold(0, 1, 1)
        with self.assertRaises(ValueError):
            q157.combined_threshold(1, 1, 1, 0, 1)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 395)
        self.assertGreater(record["closure"]["source_bindings"], 5366)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory() as td:
            target = Path(td) / "fresh.json"
            data = audit.generate()
            target.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), data)


if __name__ == "__main__":
    unittest.main()
