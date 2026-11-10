import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import q148 as q
import audit


class Tests(unittest.TestCase):
    def closure(self):
        return {
            "stage": "same-profile-compact-repair-cone-and-global-pressure-join",
            "global_identities": q.GLOBAL,
            "pressure": "Pi(X)=-integral_X^infty E^2/(2x) dx in the unchanged gauge",
            "same_exact_repair_and_inner_outer_value_cone_certified": True,
            "global_pressure_and_five_moment_identities_certified": True,
        }

    def profile(self):
        return q.finalized_profile(self.closure(), q.IDENTITY)

    def test_q147_same_profile_and_gauge_are_required(self):
        profile = self.profile()
        self.assertTrue(profile["same_profile"])
        bad = self.closure()
        bad["global_identities"] = q.GLOBAL[:-1]
        with self.assertRaises(ValueError):
            q.finalized_profile(bad, q.IDENTITY)
        bad = self.closure()
        bad["pressure"] = "pressure reset"
        with self.assertRaises(ValueError):
            q.finalized_profile(bad, q.IDENTITY)

    def test_fixed_order_linear_bound_and_supports(self):
        args = ((2, 3, 5, 7, 11), 8, 2, Q(3, 5), 7, Q(1, 1000))
        one = q.finite_order(self.profile(), 1, *args, q.IDENTITY)
        two = q.finite_order(self.profile(), 2, *args, q.IDENTITY)
        three = q.finite_order(self.profile(), 3, *args, q.IDENTITY)
        self.assertEqual(one["coefficient_norm_bound"], 3200000)
        self.assertEqual(one["support"], {"velocity": "Xvel", "pressure": "Xv", "stress": "Xb"})
        self.assertEqual(two["support"], {"velocity": "Xvel", "pressure": "Xv", "stress": "Xv"})
        self.assertEqual(set(three["support"].values()), {"Xvel"})
        self.assertTrue(one["finite_order_only"])

    def test_pressure_floor_is_sharp_and_forgets_constants(self):
        out = q.pressure_obstruction(10, 30, 3, Q(3, 7), Q(1, 5), q.IDENTITY)
        self.assertEqual(out["necessary_supremum_floor"], Q(9, 1400))
        self.assertEqual(out["independent_of"], ("K1", "K2"))
        self.assertTrue(out["sharp_after_arbitrary_X_independent_cancellation"])

    def test_fixed_order_finiteness_does_not_imply_schedule(self):
        fail = q.uniform_schedule_obligation(3, 2, 5, 7, 1, Q(1, 200), q.IDENTITY)
        self.assertFalse(fail["holds"])
        self.assertTrue(fail["not_proved_by_fixed_order_bound"])
        good = q.uniform_schedule_obligation(3, 2, 1, Q(1, 1000), Q(1, 8), Q(1, 200), q.IDENTITY)
        self.assertTrue(good["holds"])

    def test_scope_record(self):
        record = q.record()
        self.assertFalse(record["uniform_actual_pulses_complete_covariance_residual_or_blowup_certified"])
        self.assertIn("C_nm*a_n", record["uniform_blocker"])

    def test_frozen_record(self):
        with tempfile.TemporaryDirectory() as directory:
            fresh = Path(directory) / "fresh.json"
            fresh.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(fresh.read_bytes(), (audit.HERE / "evidence.json").read_bytes())


if __name__ == "__main__":
    unittest.main()
