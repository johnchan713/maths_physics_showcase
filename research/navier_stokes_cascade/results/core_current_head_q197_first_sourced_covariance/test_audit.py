import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q197


class Q197Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q197.statuses().values()))

    def test_rates_and_prefix(self):
        row = q197.rates()
        self.assertEqual(row["B_source_and_sourced_wave_log"], Q(99999, 100000))
        self.assertEqual(row["C_pressure_potential_covariance_log"], Q(149999, 100000))
        self.assertEqual(row["D_postsource_regular_log"], Q(74999, 50000))
        prefix = q197.prefix_orders()
        self.assertTrue(prefix["original_source_remains_limiting"])
        self.assertEqual(prefix["minimum"], row["B_source_and_sourced_wave_log"])

    def test_forced_pressure(self):
        row = q197.forced_pressure(7, 2, 5, 3, 2, 11)
        self.assertEqual(row["pi_over_i"], Q(5, 33))
        self.assertEqual(row["constraint_scalar"], 0)
        self.assertTrue(row["source_normal_retained"])
        self.assertTrue(row["moving_normal_retained"])

    def test_residence(self):
        self.assertEqual(q197.residence(1, 10, Q(1, 2))["residence"], 10)
        self.assertEqual(q197.residence(2, 10, Q(1, 2))["residence"], Q(2, 3))
        self.assertEqual(q197.residence(3, 10, Q(1, 2))["residence"], Q(1, 4))

    def test_degree_recurrence(self):
        row = q197.degree_recurrence((1, 2, 4, 7, 11))
        self.assertEqual(row["solution"], (2, 4, 6, 8, 12))
        self.assertTrue(row["finite_each_fixed_order"])

    def test_complete_covariance(self):
        row = q197.covariance_increment(
            {1: ((1, 2, 3), (2, 1, 4))},
            {1: ((1, 1, 1), (3, 2, 1)), 2: ((2, 3, 4), (4, 1, 2))},
        )
        self.assertTrue(row["symmetric"])
        self.assertTrue(row["all_six_present"])
        self.assertTrue(row["second_harmonic_square_retained"])
        self.assertTrue(all(row["entries"].values()))
        zero_second = q197.covariance_increment(
            {}, {1: ((1, 0, 0), (0, 0, 0)), 2: ((0, 0, 0), (0, 0, 0))},
        )
        self.assertFalse(zero_second["second_harmonic_square_retained"])

    def test_mean_pressure_and_moment(self):
        row = q197.mean_change({
            "I_DWtt_minus_DWrr_over_R": 11, "dZ_I_DWrz": 3,
            "I_R2_DWtz": 13, "I_R_DWzz": 17, "I_R_DWrr": 5,
            "I_R_DWtt": 7, "dZ_I_R2_DWrz": 19,
            "c_rho": 2, "DeltaP": 3, "dc_rho": 5, "dDeltaP": 7,
        })
        self.assertEqual(row["DeltaP_from_covariance"], Q(2197, 200))
        self.assertEqual(row["dZ_c_rho_DeltaP"], 29)
        self.assertNotEqual(row["moving_bump_term"], 0)

    def test_physical_powers_and_scope(self):
        powers = q197.physical_powers()
        self.assertEqual(powers["sourced_velocity_zero_jet"], Q(-10000001, 20000000))
        self.assertEqual(powers["postsource_residual_zero_jet"], Q(-15025001, 10000000))
        row = q197.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertFalse(row["scope"]["next_generic_signed_cycle_independently_rederived"])

    def test_invalid_inputs(self):
        for call in (q197.statuses, q197.frontier, q197.review, q197.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError): q197.rates(0)
        with self.assertRaises(ValueError): q197.forced_pressure(1, 1, 1, 0, 1, 1)
        with self.assertRaises(ValueError): q197.residence(0, 1, 1)
        with self.assertRaises(ValueError): q197.degree_recurrence(())
        with self.assertRaises(ValueError): q197.covariance_increment({0: ((1,2,3),(1,2,3))}, {})
        with self.assertRaises(ValueError): q197.mean_change({})

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 446)
        self.assertGreater(record["closure"]["edges"], 1593)
        self.assertGreater(record["closure"]["source_bindings"], 6031)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
