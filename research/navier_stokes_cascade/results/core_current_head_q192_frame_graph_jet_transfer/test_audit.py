import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q192


class Q192Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q192.statuses().values()))

    def test_inverse_recurrence(self):
        row = q192.inverse_majorants([1, 2, 3, 4, 5], Q(1, 2))
        self.assertEqual(row[0], 2)
        self.assertEqual(len(row), 5)
        self.assertTrue(all(x > 0 for x in row))

    def test_product_recurrence(self):
        left, right = [1, 2, 3, 4], [2, 3, 5, 7]
        row = q192.product_majorants(left, right)
        self.assertEqual(row[0], 2)
        self.assertEqual(row[1], 7)
        self.assertEqual(row[2], 23)

    def test_actual_acceptance_reduces_gap(self):
        r0 = Q(1, 6057072)
        for S in (64, 1000, 100000):
            length = 2 * r0 * S
            small = Q(1, 3 * S * S)
            row = q192.primitive_gap(S, length, small / 2, small / 3)
            self.assertTrue(row["closes"])
            self.assertLessEqual(row["raw_gap"], row["upper"])

    def test_auxiliary_frame_transfer(self):
        for j in range(9):
            row = q192.auxiliary_frame_transfer(j)
            self.assertEqual(row["primitive_orders_required_through"], j + 1)
            self.assertEqual(row["conservative_power"], 2 * j + 4)
            self.assertTrue(row["no_order_dependent_smallness_gate"])

    def test_corrected_graph_loss(self):
        h, kappa = Q(1, 200), Q(1, 100000)
        for a in range(8):
            for b in range(6):
                row = q192.graph_loss(a, b, h, kappa)
                self.assertEqual(row["strict_margin"], row["expected_margin"])
                self.assertTrue(row["absorbed"])
                if a + b:
                    self.assertGreater(row["strict_margin"], 0)

    def test_correction_pressure_and_frontier(self):
        fix = q192.correction()
        self.assertTrue(fix["Q191_literal_mixed_D_inequality_is_ambiguous"])
        self.assertFalse(fix["pure_polynomial_S_bound_for_physical_radial_graph_jets_claimed"])
        pressure = q192.pressure_memory()
        self.assertEqual(len(pressure["channels"]), 3)
        self.assertTrue(pressure["moving_normal_retained"])
        self.assertTrue(pressure["source_normal_retained"])
        self.assertFalse(q192.frontier()["complete_actual_field_suprema_independently_rederived_here"])

    def test_review_scope(self):
        row = q192.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["physical_graph_loss_ledger_corrected"])
        self.assertFalse(row["scope"]["complete_actual_pulse_field_bounds_rederived"])
        self.assertFalse(row["scope"]["complete_residual_improvement_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q192.statuses, q192.auxiliary_frame_transfer, q192.correction,
                     q192.pressure_memory, q192.frontier, q192.review, q192.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q192.inverse_majorants([], 1)
        with self.assertRaises(ValueError):
            q192.inverse_majorants([1], 0)
        with self.assertRaises(ValueError):
            q192.product_majorants([1], [1, 2])
        with self.assertRaises(ValueError):
            q192.primitive_gap(1, 1, 1, 1)
        with self.assertRaises(ValueError):
            q192.graph_loss(True, 0)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 442)
        self.assertGreater(record["closure"]["edges"], 1570)
        self.assertGreater(record["closure"]["source_bindings"], 5980)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
