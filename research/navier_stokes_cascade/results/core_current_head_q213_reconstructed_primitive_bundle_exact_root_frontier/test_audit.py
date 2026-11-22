import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q213


class Q213Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q213.statuses().values()))

    def test_finite_bundle(self):
        row = q213.finite_bundle([1, 9, Q(7, 2)], [Q(1, 5), Q(1, 9)], [3, 11])
        self.assertEqual(row["A_R"], 9)
        self.assertEqual(row["g_R"], Q(1, 9))
        self.assertEqual(row["K_R"], 11)

    def test_invalid_bundles(self):
        bad = [([], [1], [1]), ([1], [], [1]), ([1], [1], []),
               ([-1], [1], [1]), ([1], [0], [1]), ([True], [1], [1])]
        for args in bad:
            with self.assertRaises(ValueError):
                q213.finite_bundle(*args)

    def test_reconstructed_bundle(self):
        row = q213.reconstructed_bundle()
        self.assertTrue(row["Q212_raw_majorant_gap_inverse_obligation_closed_conditionally"])
        self.assertTrue(row["one_common_tail_and_same_profile"])
        self.assertFalse(row["numerical_high_jet_bundle_evaluated"])
        self.assertEqual(len(row["complete_pressure_memory"]), 5)

    def test_exact_frontier(self):
        row = q213.exact_frontier()
        self.assertTrue(row["named_analytical_providers_exist"])
        self.assertFalse(row["independently_rederived_in_Q213"])
        self.assertFalse(row["tests_or_hashes_are_continuum_proof"])
        self.assertIn("q_core<1", row["precise_joint_gate"])

    def test_review_scope(self):
        row = q213.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["raw_majorant_gap_inverse_bundle_joined_conditionally"])
        self.assertFalse(row["scope"]["exact_core_and_root_hypotheses_independently_reproved"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_identity_guards(self):
        for call in (q213.statuses, q213.reconstructed_bundle, q213.exact_frontier,
                     q213.review, q213.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q213.finite_bundle([1], [1], [1], identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 460)
        self.assertGreaterEqual(record["closure"]["edges"], 1700)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6200)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
