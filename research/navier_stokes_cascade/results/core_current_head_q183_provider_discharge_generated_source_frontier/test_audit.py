import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q183


class Q183Tests(unittest.TestCase):
    def test_provider_statuses_and_reconciled_frontier(self):
        self.assertTrue(all(q183.statuses().values()))
        self.assertTrue(all(q183.frontier_reconciliation().values()))

    def test_bell_recurrence_matches_exact_values(self):
        row = q183.bell_majorant((1, 2, 3), (0, 3, 5))
        self.assertEqual(row["inverse_gaussian"], [Q(1), Q(3), Q(14)])
        self.assertEqual(row["quotient"], [Q(1), Q(5), Q(29)])
        self.assertFalse(row["conditional_only"])

    def test_missing_inputs_leave_exact_gate_open(self):
        row = q183.physical_generated_source_gate()
        self.assertFalse(row["bell"]["input_suprema_supplied"])
        self.assertFalse(row["arrays_cover_requested_jet"])
        self.assertFalse(row["endpoint_flatness_proved"])
        self.assertFalse(row["gate_discharged"])
        self.assertIn("f_(j,m)/P", row["target_inequality"])

    def test_finite_inputs_are_a_conditional_control(self):
        row = q183.physical_generated_source_gate(
            raw_source_jets=(1, 2, 3),
            log_gaussian_jets=(0, 3, 5),
            endpoint_flatness_proved=True,
        )
        self.assertTrue(row["arrays_cover_requested_jet"])
        self.assertTrue(row["gate_discharged"])

    def test_review_distinguishes_route_from_proof(self):
        row = q183.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["correction"]["Q161_Bell_arithmetic_is_valid"])
        self.assertTrue(row["correction"]["Q161_synthetic_source_quotient_fixture_is_not_a_physical_supremum_proof"])
        self.assertFalse(row["scope"]["actual_physical_generated_source_supremum_proved"])
        self.assertFalse(row["scope"]["unconditional_pulse_force_exclusion_or_blowup"])

    def test_invalid_arguments(self):
        with self.assertRaises(ValueError):
            q183.statuses(identity="wrong")
        with self.assertRaises(ValueError):
            q183.bell_majorant((1,), (0, 1))
        with self.assertRaises(ValueError):
            q183.bell_majorant((1,), (-1,))
        with self.assertRaises(ValueError):
            q183.physical_generated_source_gate(harmonic=0)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 429)
        self.assertGreater(record["closure"]["edges"], 1480)
        self.assertGreater(record["closure"]["source_bindings"], 5770)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
