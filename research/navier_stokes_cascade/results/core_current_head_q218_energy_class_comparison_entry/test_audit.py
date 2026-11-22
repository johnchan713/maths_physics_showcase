import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q218


class Q218Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q218.statuses().values()))

    def test_energy_class(self):
        row = q218.energy_class()
        self.assertTrue(row["velocity_in_Linfinity_t_L2_x"])
        self.assertTrue(row["gradient_in_L2_t_L2_x"])
        self.assertTrue(row["energy_inequality"])
        self.assertFalse(row["pressure_decay_assumed"])
        self.assertFalse(row["terminal_smoothness_assumed"])

    def test_energy_entry(self):
        row = q218.energy_entry()
        self.assertEqual(row["young_inequality"], "2 <f,v> <= ||f||_2^2+E")
        self.assertTrue(row["compact_smooth_force_has_finite_rhs"])
        self.assertTrue(row["supremum_L2_required_by_Q169"])
        self.assertFalse(row["derived_from_bare_local_smoothness"])

    def test_pressure_entry(self):
        row = q218.pressure_entry()
        self.assertTrue(row["quadratic_difference_in_Linfinity_t_L1_x"])
        self.assertTrue(row["riesz_pressure_in_Hminus2"])
        self.assertTrue(row["harmonic_Hminus3_gradient_vanishes"])
        self.assertTrue(row["common_force_cancels_exactly_once"])
        self.assertFalse(row["competitor_pressure_decay_needed"])

    def test_harmonic_counterexample(self):
        row = q218.harmonic_counterexample()
        self.assertEqual(row["time_derivative"], "e_1")
        self.assertEqual(row["pressure_gradient"], "-e_1")
        self.assertEqual(row["residual"], "0")
        self.assertTrue(row["zero_initial_velocity"])
        self.assertFalse(row["finite_energy_at_positive_time"])

    def test_protected_limit(self):
        row = q218.protected_limit()
        self.assertEqual(row["velocity_exponent"], Q(101, 200))
        self.assertEqual(row["support_hole_margin"], Q(3, 101))
        self.assertTrue(row["excludes_smooth_through_terminal_energy_class_competitor"])
        self.assertFalse(row["excludes_every_locally_smooth_competitor"])

    def test_review(self):
        row = q218.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["finite_energy_comparison_premise_discharged_for_energy_class"])
        self.assertFalse(row["scope"]["bare_local_smooth_competitor_exclusion"])
        self.assertFalse(row["scope"]["unconditional_force_residual_or_blowup"])

    def test_identity_guards(self):
        for call in (q218.statuses, q218.energy_class, q218.energy_entry,
                     q218.pressure_entry, q218.harmonic_counterexample,
                     q218.protected_limit, q218.review, q218.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 468)
        self.assertGreaterEqual(record["closure"]["edges"], 1760)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6290)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
