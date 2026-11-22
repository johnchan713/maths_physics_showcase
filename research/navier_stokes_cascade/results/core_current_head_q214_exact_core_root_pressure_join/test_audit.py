import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q214


class Q214Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q214.statuses().values()))

    def test_core_gate(self):
        row = q214.core_gate(Q(2))
        self.assertEqual(row["q_core_upper"], Q(1, 10**50 * 2**127))
        self.assertTrue(row["strict_contraction"])
        self.assertTrue(row["complete_pressure_source_retained"])
        self.assertTrue(row["original_core_equations_recovered"])
        with self.assertRaises(ValueError):
            q214.core_gate(1)

    def test_root_gates(self):
        row = q214.root_gates(Q(1, 100))
        self.assertEqual(row["inverse_floors"]["compact_repair"], Q(1, 1600))
        self.assertEqual(row["inverse_floors"]["annular_five_row"], Q(1, 1220))
        self.assertEqual(row["inverse_floors"]["axial_end_correction"], Q(1, 2000))
        self.assertEqual(row["inverse_floors"]["heat_compensation"], Q(1, 100))
        self.assertTrue(row["all_displayed_floors_positive"])
        self.assertTrue(row["same_accepted_roots"])
        self.assertTrue(row["signed_variance_and_moving_phase_derivatives_positive"])
        with self.assertRaises(ValueError):
            q214.root_gates(0)

    def test_pressure_scope(self):
        row = q214.pressure_gate()
        self.assertTrue(row["Pi_core_equals_completed_leading_pressure"])
        self.assertTrue(row["smooth_endpoint_kills_homogeneous_mode"])
        self.assertEqual(len(row["positive_order_pressure_memory_retained"]), 4)
        self.assertFalse(row["full_corrected_pressure_identified_with_centrifugal_pressure"])

    def test_review(self):
        row = q214.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q213_three_part_frontier_closed_relative_to_named_providers"])
        self.assertFalse(row["scope"]["high_jet_constants_numerically_evaluated"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_identity_guards(self):
        for call in (q214.statuses, q214.core_gate, q214.root_gates,
                     q214.pressure_gate, q214.review, q214.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 464)
        self.assertGreaterEqual(record["closure"]["edges"], 1720)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6230)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
