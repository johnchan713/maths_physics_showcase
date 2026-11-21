import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q193


class Q193Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q193.statuses().values()))

    def test_sharp_Duhamel_recurrence(self):
        for n in range(16):
            self.assertEqual(q193.duhamel_degrees(n), [2 * j for j in range(n + 1)])

    def test_conservative_Duhamel_recurrence(self):
        for n in range(16):
            powers = [0] + [2 * i + 4 for i in range(1, n + 1)]
            self.assertEqual(q193.duhamel_degrees(n, powers), [7 * j for j in range(n + 1)])
            row = q193.conservative_amplitude(n)
            self.assertTrue(row["finite_polynomial_loss_for_every_fixed_order"])
            self.assertFalse(row["sharp_two_j_claimed_as_independently_needed_here"])

    def test_sourced_recurrence_is_explicitly_conditional(self):
        row = q193.sourced_degrees(3, [0, 1, 2, 3], [2, 5, 8, 11])
        self.assertEqual(row["forcing"], [2, 5, 8, 11])
        self.assertEqual(row["solution"], [3, 6, 9, 12])

    def test_all_phase_partitions(self):
        for a in range(8):
            for b in range(7):
                rows = q193.all_phase_partitions(a, b)
                self.assertEqual(len(rows), max(1, a + b))
                for row in rows:
                    self.assertTrue(row["absorbed"])
                    self.assertEqual(row["margin"], row["expected_margin"])
                    self.assertGreaterEqual(row["margin"], 0)

    def test_worst_phase_partition_saturates(self):
        row = q193.phase_partition(5, 7, 12)
        self.assertEqual(row["margin"], 0)
        self.assertEqual(row["raw_loss"], row["target_loss"])

    def test_uncut_primary_requests(self):
        for a in range(6):
            for b in range(5):
                row = q193.uncut_primary_request(a, b)
                self.assertEqual(row["amplitude_polynomial_degree"], 7 * (a + b))
                self.assertTrue(row["every_phase_partition_absorbed"])
                self.assertTrue(row["no_new_Q_loss_from_fixed_v_Duhamel"])
                self.assertTrue(row["positive_lower_bound_only_for_uncut_radial_amplitude"])

    def test_frontier_and_scope(self):
        row = q193.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertFalse(row["frontier"]["complete_cutoff_curl_pressure_potential_assembly_rederived_here"])
        self.assertFalse(row["scope"]["complete_primary_field_suprema_rederived"])
        self.assertFalse(row["scope"]["complete_residual_improvement_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q193.statuses, q193.correction, q193.frontier, q193.review, q193.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q193.duhamel_degrees(True)
        with self.assertRaises(ValueError):
            q193.duhamel_degrees(2, [0, 1])
        with self.assertRaises(ValueError):
            q193.sourced_degrees(1, [0, 1], [0])
        with self.assertRaises(ValueError):
            q193.phase_partition(1, 0, 0)
        with self.assertRaises(ValueError):
            q193.phase_partition(0, 0, 1)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 443)
        self.assertGreater(record["closure"]["edges"], 1580)
        self.assertGreater(record["closure"]["source_bindings"], 6000)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
