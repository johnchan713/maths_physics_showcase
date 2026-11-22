import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q215


class Q215Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q215.statuses().values()))

    def test_complete_request(self):
        row = q215.concrete_request()
        self.assertEqual(row["twelve_pole"], 2_626_513_221_816)
        self.assertEqual(row["safe_request"], row["twelve_pole"])
        self.assertTrue(row["derivative_ceiling_alone_is_incomplete"])
        self.assertTrue(row["finite"])

    def test_regular_exponents(self):
        self.assertEqual(q215.regular_exponent(601), Q(-401, 4_000_000))
        self.assertEqual(q215.regular_exponent(602), Q(5999, 2_500_000))
        with self.assertRaises(ValueError):
            q215.regular_exponent(-1)
        with self.assertRaises(ValueError):
            q215.regular_exponent(0, reserve=0)

    def test_bundle_scope(self):
        row = q215.bundle_join()
        self.assertTrue(row["Q207_primitive_provider_flag_closed_conditionally"])
        self.assertTrue(row["K_R_finite_including_6q_12q"])
        self.assertTrue(row["polynomial_S_degree_finite"])
        self.assertFalse(row["degree_numerically_evaluated"])

    def test_residual_memory_and_limits(self):
        row = q215.residual_join()
        self.assertEqual(len(row["six_groups"]), 6)
        self.assertTrue(row["cutoff_curls_and_pressure_derivatives_retained"])
        self.assertTrue(row["conditional_fixed_jet_residual_flatness_joined"])
        self.assertFalse(row["uniform_raw_hierarchy_constant_proved"])
        self.assertFalse(row["terminal_origin_regularized"])

    def test_review(self):
        row = q215.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q207_concrete_request_discharged_relative_to_quantified_providers"])
        self.assertFalse(row["scope"]["unconditional_residual_force_exclusion_or_blowup"])

    def test_identity_guards(self):
        for call in (q215.statuses, q215.concrete_request,
                     q215.bundle_join, q215.residual_join, q215.review, q215.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q215.regular_exponent(602, identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 465)
        self.assertGreaterEqual(record["closure"]["edges"], 1730)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6240)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
