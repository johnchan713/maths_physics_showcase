import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q189


class Q189Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q189.statuses().values()))

    def test_independent_exponent_ledger(self):
        row = q189.independent_exponent_ledger()
        self.assertTrue(all(row["checks"].values()))
        self.assertEqual(row["maximum"], 1371)
        self.assertEqual(row["envelope_exponent"], 4096)
        self.assertTrue(row["strictly_inside_envelope"])

    def test_exact_outer_ramp_bound(self):
        row = q189.rational_outer_ramp()
        self.assertEqual(row["exp_minus_three_fifths_lower"], Q(68, 125))
        self.assertEqual(row["e_minus_one_upper"], 2)
        self.assertEqual(row["normalized_p1_lower"], Q(17701, 30625))
        self.assertGreater(row["normalized_p1_lower"], Q(1, 2))
        self.assertTrue(row["above_half"])

    def test_independent_gap_join(self):
        for C in (10**9, 10**12):
            row = q189.independent_gap_join(C)
            self.assertTrue(row["all_gap_reductions_hold"])
            self.assertTrue(all(row["checks"].values()))
            self.assertEqual(row["common_required_gap"], "C^-1024")
            self.assertEqual(row["A_inverse"], "C^-4096")

    def test_same_profile_and_pressure_memory(self):
        row = q189.same_profile_contract()
        for key, value in row.items():
            if key != "endpoint_identity":
                self.assertTrue(value, key)
        self.assertEqual(row["endpoint_identity"], "Xrep=rhop*e^5=Xstar")

    def test_correction_and_scope(self):
        row = q189.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q188_envelope_and_gap_conjunction_closed_relative_to_named_analytical_providers"])
        self.assertFalse(row["scope"]["all_primitive_provider_theorems_independently_reproved_here"])
        self.assertFalse(row["scope"]["uniform_actual_pulse_or_complete_residual_proved"])

    def test_invalid_inputs(self):
        for call in (
            q189.statuses,
            q189.independent_exponent_ledger,
            q189.independent_gap_join,
            q189.same_profile_contract,
            q189.correction_and_frontier,
            q189.review,
            q189.record,
        ):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        for value in (True, 10**9 - 1, 1.0):
            with self.assertRaises(ValueError):
                q189.independent_gap_join(value)
        for args in ((), (-1,), ("1",), (True,)):
            with self.assertRaises(ValueError):
                q189._sum_exp(*args)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 439)
        self.assertGreater(record["closure"]["edges"], 1530)
        self.assertGreater(record["closure"]["source_bindings"], 5930)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
