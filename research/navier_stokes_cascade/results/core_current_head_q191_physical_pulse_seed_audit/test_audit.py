import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q191


class Q191Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q191.statuses().values()))

    def test_chart_seed(self):
        for eta in (Q(-1), Q(0), Q(3, 5), Q(1)):
            row = q191.chart_seed(eta=eta)
            self.assertEqual(row["two_D_plus_two_h"], 1)
            self.assertTrue(row["Q_cancels_from_implicit_equation"])
            self.assertGreaterEqual(row["jacobian"], Q(99, 100))

    def test_geometry_seed(self):
        for x in (Q(1), Q(33, 2), Q(1000, 3)):
            row = q191.geometry_seed(x, 64)
            self.assertTrue(row["carrier_floor"])
            self.assertTrue(row["carrier_ceiling"])
            self.assertTrue(row["Ls_upper_below_S"])

    def test_complete_pressure_constraint(self):
        row = q191.pressure_constraint()
        self.assertEqual(row["n_dot_t"], 0)
        self.assertEqual(row["constraint_derivative"], 0)
        self.assertTrue(row["moving_normal_nonzero"])
        self.assertTrue(row["source_normal_nonzero"])

    def test_Riccati_barrier(self):
        for scale in (1, 2, 1000):
            row = q191.riccati_seed(3, 5, Q(40 * scale, 3))
            self.assertTrue(row["ratio_at_most_one_eighth"])
            self.assertTrue(row["strict_outward_sign"])

    def test_physical_exponents(self):
        for a in range(5):
            for b in range(4):
                row = q191.physical_exponents(a, b)
                self.assertTrue(row["phase_losses_dominate_coefficients"])
                self.assertTrue(all(row["zero_order_scaling_checks"].values()))
                self.assertEqual(row["pressure_loss"] - row["velocity_loss"], row["sx"])

    def test_Gaussian_absorption(self):
        row = q191.gaussian_absorption(3, 5, 7, Q(1, 10))
        self.assertEqual(row["linear_exponent"], 15)
        self.assertEqual(row["maximum_exponent"], Q(1125, 2))
        self.assertTrue(row["finite_uniform_constant_exists"])

    def test_frontier_and_scope(self):
        row = q191.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertFalse(row["frontier"]["all_fixed_actual_mixed_frame_jets_independently_rederived_here"])
        self.assertFalse(row["scope"]["all_fixed_actual_physical_pulse_jets_rederived"])
        self.assertFalse(row["scope"]["unconditional_pulse_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q191.statuses, q191.chart_seed, q191.geometry_seed, q191.pressure_constraint,
                     q191.riccati_seed, q191.physical_exponents, q191.gaussian_absorption,
                     q191.frontier, q191.review, q191.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q191.chart_seed(eta=2)
        with self.assertRaises(ValueError):
            q191.geometry_seed(Q(1, 2), 1)
        with self.assertRaises(ValueError):
            q191.riccati_seed(1, 2, 8)
        with self.assertRaises(ValueError):
            q191.physical_exponents(True, 0)
        with self.assertRaises(ValueError):
            q191.gaussian_absorption(-1, 0, 0, 1)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 441)
        self.assertGreater(record["closure"]["edges"], 1560)
        self.assertGreater(record["closure"]["source_bindings"], 5980)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
