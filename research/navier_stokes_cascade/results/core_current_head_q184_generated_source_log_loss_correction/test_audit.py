import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q184


class Q184Tests(unittest.TestCase):
    def test_provider_statuses_and_analytic_derivation(self):
        self.assertTrue(all(q184.statuses().values()))
        self.assertTrue(all(q184.provider_derivation().values()))

    def test_exact_power_counterexample(self):
        row = q184.exact_power_counterexample(2, (16, 32, 64))
        self.assertEqual(row["epsilon_minus_B_times_epsilon_B_Sd"], (Q(16) ** 4, Q(32) ** 4, Q(64) ** 4))
        self.assertTrue(row["strictly_increasing"])
        self.assertFalse(row["all_band_supremum_finite"])

    def test_polynomial_absorption_is_exact_rational_bound(self):
        row = q184.polynomial_absorption(2)
        expected = Q(24) * Q(2, Q(1, 200) * Q(1, 100000)) ** 4
        self.assertEqual(row["uniform_upper"], expected)
        self.assertEqual(row["power_in_ell"], 4)

    def test_fractional_polynomial_degree_rounds_up(self):
        row = q184.polynomial_absorption(Q(9, 4))
        self.assertEqual(row["degree_ceiling"], 3)
        self.assertEqual(row["power_in_ell"], 6)

    def test_zero_degree_needs_no_loss_constant(self):
        row = q184.polynomial_absorption(0)
        self.assertEqual(row["uniform_upper"], 1)

    def test_corrected_rates(self):
        row = q184.corrected_rates()
        self.assertEqual(row["B_log"], Q(3, 2) - 2 * Q(1, 100000))
        self.assertEqual(row["B_pure"], Q(3, 2) - 3 * Q(1, 100000))
        self.assertEqual(row["delta_log"], Q(1, 2) - 3 * Q(1, 100000))
        self.assertEqual(row["delta_pure"], Q(1, 2) - 4 * Q(1, 100000))
        self.assertEqual(row["target_pure_margin_over_required"], Q(1, 50000))
        self.assertTrue(row["pure_source_gain_positive"])

    def test_fixed_cycle_contract_distinguishes_classes(self):
        row = q184.fixed_cycle_source_contract(polynomial_degree=3)
        self.assertTrue(row["same_domain_raw_source_Ai_exist_conditionally"])
        self.assertTrue(row["same_domain_log_P_Ki_exist_conditionally"])
        self.assertTrue(row["smooth_endpoint_extension_exists_conditionally"])
        self.assertFalse(row["Q183_exact_power_bound_follows_without_margin"])
        zero = q184.fixed_cycle_source_contract(polynomial_degree=0)
        self.assertTrue(zero["Q183_exact_power_bound_follows_without_margin"])

    def test_review_correction_and_limits(self):
        row = q184.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["correction"]["Q183_misidentified_missing_Ai_Ki_as_first_conditional_gap"])
        self.assertTrue(row["scope"]["log_renormalized_fixed_cycle_source_contract_closed_conditionally"])
        self.assertFalse(row["scope"]["original_exact_delta_pure_power_contract_proved"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            q184.statuses(identity="wrong")
        for bad in (-1, True, 0.5):
            with self.assertRaises(ValueError):
                q184.polynomial_absorption(bad)
        with self.assertRaises(ValueError):
            q184.polynomial_absorption(1, eta=0)
        with self.assertRaises(ValueError):
            q184.exact_power_counterexample(0)
        with self.assertRaises(ValueError):
            q184.corrected_rates(-1)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 430)
        self.assertGreater(record["closure"]["edges"], 1490)
        self.assertGreater(record["closure"]["source_bindings"], 5780)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
