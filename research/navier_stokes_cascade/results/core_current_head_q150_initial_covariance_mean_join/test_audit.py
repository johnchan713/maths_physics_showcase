import json
import tempfile
import unittest
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import audit
import q150 as q


class Tests(unittest.TestCase):
    def q149(self):
        path = q.HERE.parent / "core_current_head_q149_primary_pulse_join/evidence.json"
        return json.loads(path.read_text())["symbolic"]

    def test_same_physical_primary_initial_join(self):
        out = q.initial_join(self.q149(), q.IDENTITY)
        self.assertIn("all six", out["covariance"])
        self.assertIn("Pi1", out["pressure"])
        self.assertEqual(out["initial_recomputed_defects"], "P,Jtheta,Jz belong to S_2 and are generally nonzero")
        self.assertFalse(out["later_signed_dynamic_cycle_certified"])

    def test_missing_parent_or_provider_contract_rejected(self):
        bad = self.q149()
        bad["status"] = "primary-only"
        with self.assertRaises(ValueError):
            q.initial_join(bad, q.IDENTITY)
        original = q.provider_flags
        try:
            flags = deepcopy(original())
            flags["complete_covariance_initial_means"]["complete_axial_velocity_and_moving_pressure_moment_retained"] = False
            q.provider_flags = lambda: flags
            with self.assertRaises(ValueError):
                q.initial_join(self.q149(), q.IDENTITY)
        finally:
            q.provider_flags = original

    def test_initial_gain_is_positive_but_initial_only(self):
        out = q.initial_rates(Q(1, 200))
        self.assertEqual(out["epsilon_initial_residual_order"], Q(199997, 100000))
        self.assertEqual(out["later_desired_order_D"], Q(74999, 50000))
        self.assertEqual(out["initial_margin_over_D"], Q(49999, 100000))
        self.assertTrue(out["positive_initial_margin"])
        self.assertFalse(out["transfer_to_later_nonzero_mean_state"])

    def test_q_conversion_of_initial_residual_order(self):
        out = q.initial_rates(Q(1, 200))
        self.assertEqual(out["Q_initial_residual_exponent"], Q(199997, 20000000))
        self.assertGreater(out["Q_initial_residual_exponent"], 0)

    def test_precise_later_feedback_gap(self):
        out = q.frontier()
        self.assertEqual(out["B"], Q(99999, 100000))
        self.assertEqual(out["C"], Q(149999, 100000))
        self.assertEqual(out["D"], Q(74999, 50000))
        self.assertEqual(out["missing_gain"], Q(49999, 100000))
        self.assertFalse(out["coarse_mean_size_alone_suffices"])
        self.assertIn("Pperp", out["inequality"])

    def test_invalid_identity_or_ranges_rejected(self):
        with self.assertRaises(ValueError):
            q.initial_join(self.q149(), "other")
        for h in (0, Q(1, 100), -1):
            with self.assertRaises(ValueError):
                q.initial_rates(h)
        with self.assertRaises(ValueError):
            q.frontier(identity="other")

    def test_scope_record(self):
        record = q.record()
        self.assertFalse(record["later_signed_dynamic_cycle_or_complete_residual_certified"])
        self.assertFalse(record["blowup_certified"])
        self.assertIn("initial zero-old-mean stage only", record["scope"])
        self.assertIn("1/2-kappa", record["obstruction"])

    def test_frozen_record(self):
        with tempfile.TemporaryDirectory() as directory:
            fresh = Path(directory) / "fresh.json"
            fresh.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(fresh.read_bytes(), (audit.HERE / "evidence.json").read_bytes())


if __name__ == "__main__":
    unittest.main()
