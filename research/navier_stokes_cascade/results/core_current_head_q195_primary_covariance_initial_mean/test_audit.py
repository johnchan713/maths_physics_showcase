import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q195


class Q195Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q195.statuses().values()))

    def test_complete_covariance(self):
        row = q195.covariance((1, 2, 3), (4, 5, 6))
        self.assertTrue(row["symmetric"])
        self.assertTrue(row["all_six_present"])
        self.assertEqual(row["entries"]["rtheta"], Q(11))
        self.assertEqual(row["entries"]["thetaz"], Q(18))
        self.assertEqual(row["mixed_sine_cosine_mean"], 0)

    def test_physical_covariance_loss(self):
        for a in range(7):
            for b in range(6):
                v = q195.physical_loss("velocity", a, b)
                w = q195.physical_loss("covariance", a, b)
                p = q195.physical_loss("wave_pressure", a, b)
                self.assertEqual(w["Q_loss"], p["Q_loss"])
                self.assertEqual(w["Q_loss"] - v["Q_loss"], Q(201, 400))
        self.assertEqual(q195.physical_loss("covariance")["shell_weight"], "zeta")
        self.assertEqual(q195.physical_loss("covariance")["Gaussian_weight"], "P^2")

    def test_same_family_match(self):
        row = q195.principal_match(
            ((Q(1, 4), 6, 10, 2), (Q(3, 4), 2, 6, 2)), 6, 14
        )
        self.assertTrue(row["exact"])
        self.assertEqual(row["cover"], 1)
        bad = q195.principal_match(((1, 6, 10, 2),), 13, 20)
        self.assertFalse(bad["exact"])

    def test_initial_target_component_ledger(self):
        row = q195.initial_targets({
            "I_Wtt_minus_Wrr_over_R": 11,
            "dZ_I_Wrz": 3,
            "I_R2_Wtz": 13,
            "I_R_Wzz": 17,
            "I_R_Wrr": 5,
            "I_R_Wtt": 7,
            "dZ_I_R2_Wrz": 19,
        })
        self.assertEqual(row["P"], Q(2197, 200))
        self.assertEqual(row["Jtheta"], 13)
        self.assertEqual(row["Jz"], Q(4419, 400))
        self.assertEqual(set(row["used_components"]) | {"rtheta"}, set(q195.COMPONENTS))

    def test_whole_moving_pressure_product(self):
        row = q195.moving_pressure_derivative(2, 3, 5, 7)
        self.assertEqual(row["whole_product"], 29)
        self.assertEqual(row["bump_motion_term"], 15)
        self.assertTrue(row["omission_changes_identity"])

    def test_rates_and_scope(self):
        rates = q195.rates()
        self.assertEqual(rates["covariance_remainder"], Q(99999, 50000))
        self.assertEqual(rates["initial_tangential_residual"], Q(199997, 100000))
        row = q195.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["complete_real_primary_covariance_rederived"])
        self.assertFalse(row["scope"]["later_complex_cycles_and_generated_sources_rederived"])

    def test_invalid_inputs(self):
        for call in (q195.statuses, q195.frontier, q195.review, q195.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q195.covariance((1, 2), (3, 4, 5))
        with self.assertRaises(ValueError):
            q195.physical_loss("other")
        with self.assertRaises(ValueError):
            q195.principal_match((), 0, 0)
        with self.assertRaises(ValueError):
            q195.initial_targets({})

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 444)
        self.assertGreater(record["closure"]["edges"], 1580)
        self.assertGreater(record["closure"]["source_bindings"], 6000)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
