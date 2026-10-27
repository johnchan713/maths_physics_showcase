from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from corner_replay import cartesian_constant, closed_chart, residual, review, target


class LogReservedCartesianCornerReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_closed_chart_exact_geometry(self):
        chart = closed_chart()
        self.assertEqual(chart["terminal_identity"], 1)
        self.assertEqual(chart["jacobian_floor"], Q(99, 100))
        self.assertEqual(chart["normalized_support_margin"], Q(3, 8012))

    def test_scale_incompatible_floors_are_absent(self):
        chart = closed_chart()
        self.assertFalse(chart["physical_support_floor_used"])
        self.assertFalse(chart["edge_weight_floor_used"])

    def test_cartesian_factor(self):
        for power in (0, 1, 2, 9, 33):
            self.assertEqual(cartesian_constant(power), 4 ** power)

    def test_default_target_uses_log_reserved_stage(self):
        row = target()
        parent = residual.comparison()
        self.assertEqual(row["comparison_stage"], 12817)
        self.assertEqual(row["comparison_stage"], parent["log_reserved_comparison_stage"])
        self.assertGreaterEqual(row["closed_chart_exponent"], 1)

    def test_all_six_groups_and_pressure_memory(self):
        row = target()
        self.assertEqual(tuple(row["residual_groups"]), residual.RESIDUAL_GROUPS)
        self.assertIn("all six covariance entries", row["pressure_memory"])
        self.assertIn("exceptional compact tails", row["pressure_memory"])

    def test_multiple_jet_and_power_targets(self):
        for h in (Q(1, 200), Q(9, 1000), Q(1, 10000)):
            for m, power in ((0, 0), (0, 10), (3, 5), (8, 17), (20, 2)):
                row = target(m, power, h)
                self.assertTrue(all(row["checks"].values()))
                self.assertGreaterEqual(row["retained_exponent_margin"], 0)

    def test_invalid_inputs_rejected(self):
        for value in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                target(value, 0)
        with self.assertRaises(ValueError):
            closed_chart(Q(1, 100))

    def test_review_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        scope = result["scope"]
        self.assertTrue(scope["same_field_all_fixed_residual_jets_cartesian_corner_flat_conditionally"])
        self.assertFalse(scope["localization_force_energy_comparison_replayed_here"])
        self.assertFalse(scope["unconditional_pulse_residual_force_or_blowup_claimed"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["unconditional_pulse_residual_force_or_blowup_claimed"] = True
        with self.assertRaises(ValueError):
            audit.validate(bad)

    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
