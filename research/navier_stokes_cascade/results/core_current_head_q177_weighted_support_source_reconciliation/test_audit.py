import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q177


class Q177Tests(unittest.TestCase):
    def test_statuses(self):
        self.assertTrue(all(q177.statuses().values()))

    def test_frequency_acceptance_gate(self):
        self.assertFalse(q177.frequency_acceptance_gate(2)["log_threshold_below_32_G4"])
        row = q177.frequency_acceptance_gate(4)
        self.assertTrue(row["log_threshold_below_32_G4"])
        self.assertTrue(row["cone_absorption"])
        self.assertEqual(row["log_threshold_slack"], 16 * 4 ** 4 - 193 * 4)
        with self.assertRaises(ValueError):
            q177.frequency_acceptance_gate(0)

    def test_support_geometry(self):
        row = q177.support_geometry()
        self.assertTrue(row["strict_separation"])
        self.assertEqual(row["reserved_gap"], 3)
        with self.assertRaises(ValueError):
            q177.support_geometry(i4_left=0)

    def test_support_geometry_rejects_overlap(self):
        self.assertFalse(q177.support_geometry(i3_right=4)["strict_separation"])
        self.assertFalse(q177.support_geometry(i4_left=Q(5, 2))["strict_separation"])

    def test_linear_inverse_gate(self):
        row = q177.linear_inverse_gate()
        self.assertEqual(row["E_inverse_upper"], Q(393781, 500))
        self.assertTrue(row["strict_inverse_gate"])
        self.assertFalse(q177.linear_inverse_gate(E_inverse=800)["strict_inverse_gate"])
        with self.assertRaises(ValueError):
            q177.linear_inverse_gate(U_inverse=0)

    def test_outer_quotient_jets(self):
        row = q177.outer_quotient_jets(13)
        self.assertEqual(row["cubic_coefficient"], Q(13, 8))
        self.assertEqual(row["third_delta_derivative"], Q(39, 4))
        self.assertEqual(row["second_eta_derivative"], 52)
        with self.assertRaises(ValueError):
            q177.outer_quotient_jets(0)

    def test_reserved_pressure(self):
        row = q177.reserved_pressure(2, 5, 7, 11)
        self.assertEqual(row, {"Pi1": Q(34, 5), "Pi2": Q(11), "higher_pressure": Q(0)})
        with self.assertRaises(ValueError):
            q177.reserved_pressure(1, 0, 0, 0)

    def test_review_correction_and_scope(self):
        row = q177.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["correction"]["Q176_missing_selected_supremum_route_was_too_broad"])
        self.assertTrue(row["correction"]["external_independent_validation_of_exact_core_and_finite_lemmas_still_open"])
        self.assertFalse(row["correction"]["historical_frozen_evidence_modified"])
        self.assertTrue(row["established"]["finite_coefficient_and_reserved_velocity_supports_with_pressure_memory"])
        self.assertFalse(row["scope"]["unconditional_uniform_actual_pulses"])
        self.assertFalse(row["scope"]["complete_covariance_mean_residual_or_blowup"])

    def test_invalid_identity(self):
        for fn in (q177.statuses, q177.review):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 422)
        self.assertGreater(record["closure"]["source_bindings"], 5690)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
