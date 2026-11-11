import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q179


class Q179Tests(unittest.TestCase):
    def test_statuses(self):
        self.assertTrue(all(q179.statuses().values()))

    def test_exact_rates(self):
        row = q179.rates()
        self.assertEqual(row["B"], Q(99999, 100000))
        self.assertEqual(row["C"], Q(149999, 100000))
        self.assertEqual(row["D"], Q(74999, 50000))
        self.assertEqual(row["initial_residual_order"], Q(199997, 100000))
        self.assertEqual(row["Q_initial_residual_exponent"], Q(199997, 20000000))
        self.assertEqual(row["initial_margin_over_D"], Q(49999, 100000))
        self.assertEqual(row["complete_curl_margin"], Q(1, 2))
        self.assertEqual(row["first_feedback_order"], row["D"])
        with self.assertRaises(ValueError):
            q179.rates(kappa=0)

    def test_principal_assembly(self):
        good = [(Q(1, 4), 2, 3, 5, 15), (Q(3, 4), 2, 3, 5, 15)]
        self.assertTrue(q179.principal_assembly(good, 30)["exact"])
        self.assertFalse(q179.principal_assembly(good[:1], 30)["exact"])
        with self.assertRaises(ValueError):
            q179.principal_assembly([], 0)

    def test_two_moment_rows(self):
        self.assertTrue(q179.moment_gate(7, 7, 11, 11)["both_exact"])
        self.assertFalse(q179.moment_gate(7, 6, 11, 11)["both_exact"])
        self.assertFalse(q179.moment_gate(7, 7, 11, 10)["both_exact"])

    def test_review_correction_and_scope(self):
        row = q179.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["established"]["initial_complete_six_entry_covariance_mean_cancellation"])
        self.assertTrue(row["established"]["same_physical_family_principal_stress_and_first_feedback_D"])
        self.assertTrue(row["correction"]["Q178_raw_negative_covariance_power_is_not_the_initial_stage_frontier"])
        self.assertFalse(row["obstruction"]["ordinary_unweighted_compact_seminorms_suffice"])
        self.assertFalse(row["scope"]["P3_universal_finite_cycle_target_gate_discharged"])
        self.assertFalse(row["scope"]["complete_residual_force_exclusion_or_blowup"])

    def test_invalid_identity(self):
        for fn in (q179.statuses, q179.review):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 424)
        self.assertGreater(record["closure"]["source_bindings"], 5720)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
