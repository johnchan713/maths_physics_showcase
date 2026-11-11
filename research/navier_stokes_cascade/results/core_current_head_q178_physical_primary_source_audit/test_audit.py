import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q178


class Q178Tests(unittest.TestCase):
    def test_statuses(self):
        self.assertTrue(all(q178.statuses().values()))

    def test_physical_losses(self):
        zero = q178.physical_losses(0, 0)
        self.assertEqual(zero["spatial_loss"], Q(201, 400))
        self.assertEqual(zero["time_loss"], Q(403, 400))
        self.assertEqual(zero["velocity_growth"], Q(201, 400))
        self.assertEqual(zero["pressure_growth"], Q(201, 200))
        mixed = q178.physical_losses(2, 1)
        self.assertEqual(mixed["potential_growth"], Q(161, 80))
        with self.assertRaises(ValueError):
            q178.physical_losses(-1, 0)

    def test_borel_tail(self):
        row = q178.borel_tail(5, Q(3, 4))
        self.assertEqual(row["exact"], (Q(3, 8) ** 5) / Q(5, 8))
        self.assertTrue(row["bounded"])
        self.assertEqual(q178.borel_tail(1, 1)["uniform_upper"], 1)
        with self.assertRaises(ValueError):
            q178.borel_tail(0, 1)

    def test_covariance_growth(self):
        early = q178.covariance_growth(10)
        later = q178.covariance_growth(20)
        self.assertEqual(early["Q_exponent"], Q(-10000001, 10000000))
        self.assertTrue(early["grows_with_band"])
        self.assertFalse(early["absolutely_decaying"])
        self.assertGreater(later["binary_growth_exponent"], early["binary_growth_exponent"])
        with self.assertRaises(ValueError):
            q178.covariance_growth(-1)

    def test_review_boundary_and_scope(self):
        row = q178.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["established"]["conditional_uniform_actual_primary_pulse_values_fixed_derivatives_and_supports"])
        self.assertTrue(row["correction"]["uncut_radial_amplitude_lower_bound_is_not_a_pointwise_cutoff_velocity_lower_bound"])
        self.assertTrue(row["correction"]["finite_primary_bounds_do_not_verify_generated_source_membership"])
        self.assertFalse(row["scope"]["all_generated_nonlinear_source_inputs_verified"])
        self.assertFalse(row["scope"]["complete_covariance_mean_residual_or_blowup"])

    def test_invalid_identity(self):
        for fn in (q178.statuses, q178.review):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 423)
        self.assertGreater(record["closure"]["source_bindings"], 5700)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
