import json
import tempfile
import unittest
from pathlib import Path
import audit
import q221


class Q221Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q221.statuses().values()))

    def test_exact_integer_norms(self):
        for k1 in range(-12, 13):
            for k2 in range(-12, 13):
                if (k1, k2) == (0, 0):
                    continue
                row = q221.quadratic_norms(k1, k2)
                self.assertTrue(row["both_norms_nonzero"])
                self.assertGreaterEqual(abs(row["vr_integer_norm"]), 1)
                self.assertGreaterEqual(abs(row["vt_integer_norm"]), 1)
        with self.assertRaises(ValueError):
            q221.quadratic_norms(0, 0)

    def test_fourier_contract(self):
        for order in (0, 1, 17, 2_400_066):
            row = q221.fourier_contract(order)
            self.assertEqual(row["source_order"], order + 4)
            self.assertEqual(row["post_integration_per_mode_power"], -3)
            self.assertTrue(row["zero_Haar_required"])
            self.assertFalse(row["giant_constant_materialized"])

    def test_fast_contract(self):
        row = q221.fast_contract(17)
        self.assertEqual(row["tangential_axial_source_order"], 21)
        self.assertEqual(row["radial_potential_source_order"], 22)
        self.assertEqual(row["whole_axial_pressure_source_order"], 26)
        self.assertTrue(row["actual_A1_axial_repair_retained"])

    def test_slow_contract(self):
        pole = 2_626_513_221_816
        row = q221.slow_contract(2_400_066, pole)
        self.assertEqual(row["output_pole"], pole)
        self.assertTrue(row["five_moment_rows_exact"])
        self.assertTrue(row["whole_convex_potential_hull_used"])
        self.assertTrue(row["interior_weight_floor_preserves_every_fixed_pole"])

    def test_axial_conditioning_obstruction(self):
        row = q221.axial_conditioning()
        self.assertTrue(row["finite_for_selected_lambda"])
        self.assertFalse(row["uniform_as_lambda_downarrow_zero"])
        self.assertFalse(row["angular_block_singular_at_lambda_zero"])

    def test_pressure_memory(self):
        row = q221.pressure_memory()
        self.assertTrue(row["c_rho_Z_P_retained"])
        self.assertTrue(row["actual_axial_velocity_retained"])
        self.assertTrue(row["all_six_covariances_retained_in_inputs"])
        self.assertFalse(row["operator_creates_its_own_source_bounds"])

    def test_review(self):
        row = q221.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q216_remaining_inverse_quantifiers_supported"])
        self.assertFalse(row["scope"]["lambda_uniform_slow_inverse_claim"])
        self.assertFalse(row["scope"]["unconditional_force_exclusion_or_blowup"])

    def test_identity_and_type_guards(self):
        for call in (q221.statuses, q221.axial_conditioning,
                     q221.pressure_memory, q221.review, q221.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        for call, args in ((q221.quadratic_norms, (True, 1)),
                           (q221.fourier_contract, (-1,)),
                           (q221.fast_contract, (1.5,)),
                           (q221.slow_contract, (1, -1))):
            with self.assertRaises(ValueError):
                call(*args)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 471)
        self.assertGreaterEqual(record["closure"]["edges"], 1788)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6331)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
