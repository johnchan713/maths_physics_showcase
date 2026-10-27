from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from endpoint_compatibility import corner, endpoint, localization, request, review


class LogReservedEndpointCompatibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_default_order_join(self):
        row = request()
        self.assertEqual(row["punctured_field_orders"], {
            "potential": 4, "direct_swirl": 3, "pressure": 2,
        })
        self.assertEqual(row["localization_Cauchy_orders"], row["punctured_field_orders"])
        self.assertEqual(row["residual_Cauchy_jet"], 1)

    def test_extra_time_derivative_raises_only_finite_order(self):
        for a, j in ((0, 0), (1, 3), (20, 4), (6, 70), (100, 100)):
            row = request(a, j)
            total = a + j
            self.assertEqual(row["residual_Cauchy_jet"], total + 1)
            self.assertEqual(row["punctured_field_orders"]["potential"], total + 4)
            self.assertTrue(all(row["checks"].values()))

    def test_request_stage_estimates_same_field(self):
        low, high = request(0, 0, 0), request(20, 20, 20)
        self.assertGreater(high["comparison_stage"], low["comparison_stage"])
        self.assertTrue(low["checks"]["same_fixed_field_is_used"])
        self.assertTrue(high["checks"]["same_fixed_field_is_used"])

    def test_exact_parent_ledgers_match(self):
        for m in (0, 1, 5, 100):
            force = localization.force_requirements(m)
            punctured = endpoint.requirements(m)
            self.assertEqual(force["Cauchy_potential_order"], punctured["potential_for_future_residual_Cauchy"])
            self.assertEqual(force["Cauchy_direct_swirl_order"], punctured["direct_for_future_residual_Cauchy"])
            self.assertEqual(force["Cauchy_pressure_order"], punctured["pressure_for_future_residual_Cauchy"])

    def test_corner_parent_retains_six_groups(self):
        row = corner.target(1, 3)
        self.assertEqual(tuple(row["residual_groups"]), corner.residual.RESIDUAL_GROUPS)
        self.assertGreaterEqual(row["retained_exponent_margin"], 0)

    def test_invalid_inputs_rejected(self):
        for value in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                request(value, 0, 0)

    def test_review_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        scope = result["scope"]
        self.assertTrue(scope["closed_time_localization_inputs_joined_conditionally"])
        self.assertFalse(scope["terminal_velocity_or_pressure_regular_at_origin_proved"])
        self.assertFalse(scope["localization_and_compact_force_replayed_here"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["terminal_velocity_or_pressure_regular_at_origin_proved"] = True
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
