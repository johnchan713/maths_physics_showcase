import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q223


class Q223Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q223.statuses().values()))

    def test_coefficient_target(self):
        self.assertEqual(q223.coefficient_target(3, 5, 7, 11, 13, Q(2, 9)),
                         (Q(3), Q(-18), Q(5), Q(-11), Q(13)))
        with self.assertRaises(ValueError):
            q223.coefficient_target(1, 2, 3, 4, 5, 0)

    def test_linear_contract(self):
        row = q223.coefficient_contract()
        self.assertFalse(row["linear_solve_requires_small_target"])
        self.assertFalse(row["uniform_in_order_or_derivative"])

    def test_support_ledger(self):
        self.assertEqual(q223.support_ledger(1),
                         {"velocity": "Xvel", "pressure": "Xv", "stress": "Xb"})
        self.assertEqual(q223.support_ledger(2),
                         {"velocity": "Xvel", "pressure": "Xv", "stress": "Xv"})
        self.assertEqual(q223.support_ledger(9),
                         {"velocity": "Xvel", "pressure": "Xvel", "stress": "Xvel"})
        with self.assertRaises(ValueError):
            q223.support_ledger(0)

    def test_reserved_geometry(self):
        row = q223.reserved_geometry()
        self.assertTrue(row["strict_separation"])
        self.assertTrue(row["positive_order_velocity_zero_on_reserved_interval"])
        self.assertFalse(row["pressure_inferred_zero_from_velocity_support"])

    def test_minimax_constant_floor_is_sharp(self):
        row = q223.minimax_constant_floor(Q(-9, 7), Q(5, 11))
        self.assertEqual(abs(row["left_at_optimum"]), row["sharp_floor"])
        self.assertEqual(abs(row["right_at_optimum"]), row["sharp_floor"])
        self.assertEqual(row["left_at_optimum"], -row["right_at_optimum"])

    def test_pressure_two_point_floor(self):
        row = q223.pressure_two_point_floor(3, 5, 7, 11)
        expected = Q(3 * 5, 8) * (Q(1, 7) - Q(1, 11))
        self.assertEqual(row["sharp_floor"], expected)
        self.assertTrue(row["strictly_positive"])
        with self.assertRaises(ValueError):
            q223.pressure_two_point_floor(3, 5, 11, 7)

    def test_Pi2_negative_control(self):
        row = q223.minimax_constant_floor(Q(13, 17), Q(13, 17))
        self.assertEqual(row["sharp_floor"], 0)
        claim = q223.corrected_claim()
        self.assertFalse(claim["joint_Pi1_Pi2_positive_floor"])
        self.assertFalse(claim["Pi2_deletion_permitted"])

    def test_frontier(self):
        row = q223.frontier()
        self.assertEqual(row["rule_selected_Borel_schedule"], "already closed by Q204")
        self.assertFalse(row["uniform_raw_uncut_generated_hierarchy_constant"])
        self.assertFalse(row["complete_covariance_pressure_viscosity_cutoff_interaction_residual"])

    def test_review(self):
        row = q223.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["sharp_Pi1_radial_pressure_floor_derived"])
        self.assertFalse(row["scope"]["positive_Pi2_floor_claimed"])
        self.assertFalse(row["scope"]["complete_residual_or_blowup"])

    def test_identity_guards(self):
        calls = (
            q223.statuses, q223.coefficient_contract, q223.reserved_geometry,
            q223.reserved_pressure_ledger, q223.corrected_claim, q223.frontier,
            q223.review, q223.record,
        )
        for call in calls:
            with self.assertRaises(ValueError):
                call(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 473)
        self.assertGreaterEqual(record["closure"]["edges"], 1800)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6350)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
