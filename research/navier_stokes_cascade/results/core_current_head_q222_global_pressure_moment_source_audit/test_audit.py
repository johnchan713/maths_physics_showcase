import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q222


class Q222Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q222.statuses().values()))

    def test_finite_transforms_round_trip(self):
        row = tuple(map(Q, (3, 5, 7, 11, 13)))
        for eta in (Q(0), Q(2, 9), Q(-7, 13)):
            self.assertEqual(q222.annular_inverse(q222.annular_transform(row, eta), eta), row)
        for lam in (Q(1, 100), Q(5, 17), Q(19)):
            self.assertEqual(q222.compact_inverse(q222.compact_transform(row, lam), lam), row)
        with self.assertRaises(ValueError):
            q222.compact_transform(row, 0)

    def test_global_heat_common_tail(self):
        self.assertEqual(q222.global_heat_join(0, 0)["delta_R"], 0)
        self.assertNotEqual(q222.global_heat_join(0, Q(1, 19))["delta_R"], 0)
        self.assertFalse(q222.global_heat_join(0, 0)["endpoint_I_is_global_R"])

    def test_leading_global_ledger(self):
        row = q222.leading_global_ledger()
        self.assertFalse(row["intermediate_pressure_reset"])
        self.assertTrue(row["all_five_rows_required_for_exterior_stress_join"])
        self.assertTrue(row["pressure_derivative_defect_is_independent"])

    def test_compact_pressure_moment(self):
        row = q222.compact_pressure_moment(7, 11, Q(13, 17))
        self.assertEqual(row["integral_R_pbar"], -Q(7, 2) + Q(13, 17) * 11)
        self.assertFalse(row["A0_pointwise_remainder_deleted"])

    def test_moving_product_derivative(self):
        self.assertEqual(q222.moving_product_derivative(2, 3, 5, 7, 11, 13),
                         13 * (2 + 3 * 5 + 7 * 11))
        self.assertNotEqual(q222.moving_product_derivative(2, 0, 5, 7, 11, 13),
                            q222.moving_product_derivative(2, 3, 5, 7, 11, 13))

    def test_mean_moment_ledger(self):
        row = q222.mean_moment_ledger()
        self.assertTrue(row["whole_moving_product_differentiated"])
        self.assertTrue(row["compact_radial_boundary_terms_zero"])
        self.assertFalse(row["pointwise_rows_determined_by_moments"])

    def test_preservation_ledger(self):
        row = q222.preservation_ledger()
        self.assertFalse(row["nonzero_harmonic_waves_change_Haar_velocity_moments"])
        self.assertFalse(row["zero_Haar_fast_curl_change_Haar_velocity_moments"])
        self.assertFalse(row["zero_Haar_correlations_discarded"])
        self.assertTrue(row["Pi1_Pi2_cancel_only_in_same_complete_base_subtraction"])

    def test_projected_targets(self):
        row = q222.projected_targets()
        self.assertTrue(row["h1_zero_moment"])
        self.assertTrue(row["h2_zero_moment"])
        self.assertFalse(row["projection_deletes_pressure_or_mean_scalar"])
        self.assertTrue(row["two_edge_primitive_right_collar_requires_exact_zero_moment"])

    def test_review(self):
        row = q222.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["leading_five_global_rows_joined_on_common_tail"])
        self.assertFalse(row["scope"]["pointwise_Haar_rows_inferred_from_moments"])
        self.assertFalse(row["scope"]["unconditional_force_exclusion_or_blowup"])

    def test_identity_guards(self):
        for call in (q222.statuses, q222.leading_global_ledger,
                     q222.mean_moment_ledger, q222.preservation_ledger,
                     q222.projected_targets, q222.review, q222.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 472)
        self.assertGreaterEqual(record["closure"]["edges"], 1798)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6345)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
