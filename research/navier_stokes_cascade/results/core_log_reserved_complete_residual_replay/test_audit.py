from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from replay import (
    RESIDUAL_GROUPS, actual_source, compact_reserve_accounts, comparison,
    diagonal, logarithmic_absorption, review,
)


class LogReservedCompleteResidualReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_zero_reserve_counterexample_remains_unbounded(self):
        self.assertEqual(actual_source.no_log_counterexample(17)["pure_power_ratio"], 18)
        self.assertEqual(actual_source.no_log_counterexample(1700)["pure_power_ratio"], 1701)

    def test_absorption_requires_positive_reserve(self):
        self.assertEqual(logarithmic_absorption(2, Q(1, 5)), 5)
        for reserve in (0, -1):
            with self.assertRaises(ValueError):
                logarithmic_absorption(2, reserve)

    def test_compact_accounts_are_separate(self):
        accounts = compact_reserve_accounts()
        expected = Q(1, 200) * diagonal.KAPPA / 2
        self.assertIsNot(accounts["A0"], accounts["A1"])
        self.assertEqual(accounts["A0"]["retained_strict_margin"], expected)
        self.assertEqual(accounts["A1"]["spent_on_own_source_log"], expected)
        self.assertFalse(accounts["A0"]["may_pay_regular_or_cutoff_log"])

    def test_one_stage_supplies_the_regular_log_reserve(self):
        row = comparison()
        self.assertEqual(row["log_reserved_comparison_stage"],
                         row["base_comparison_stage"] + 1)
        self.assertEqual(row["raw_regular_exponent"] - row["base_regular_exponent"],
                         row["regular_stage_gain"])
        self.assertEqual(row["regular_log_spend"], row["regular_retained_margin"])

    def test_default_stage_and_all_checks(self):
        row = comparison()
        self.assertEqual(row["base_comparison_stage"], 12816)
        self.assertEqual(row["log_reserved_comparison_stage"], 12817)
        self.assertTrue(all(row["checks"].values()))

    def test_six_exact_residual_groups(self):
        row = comparison()
        self.assertEqual(tuple(row["residual_group_exponents"]), RESIDUAL_GROUPS)
        self.assertEqual(row["residual_group_exponents"]["old_new"],
                         row["residual_group_exponents"]["new_old"])
        self.assertIn("Delta e", row["identity"])
        self.assertIn("grad pi", row["identity"])

    def test_multiple_jet_and_power_targets(self):
        for h in (Q(1, 200), Q(9, 1000), Q(1, 10000)):
            for m, target in ((0, 0), (1, 3), (5, 11), (20, 2)):
                row = comparison(m, target, h)
                self.assertTrue(all(row["checks"].values()))
                self.assertGreaterEqual(row["complete_minimum_exponent"], target + 1)

    def test_invalid_inputs_rejected(self):
        for value in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                comparison(value, 0)
        with self.assertRaises(ValueError):
            compact_reserve_accounts(kappa=0)

    def test_review_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        scope = result["scope"]
        self.assertTrue(scope["log_reserved_six_group_complete_local_residual_flat_conditionally"])
        self.assertFalse(scope["unconditional_pulse_force_comparison_or_blowup_claimed"])
        self.assertFalse(scope["localization_force_energy_and_comparison_replayed_after_source_correction"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["unconditional_pulse_force_comparison_or_blowup_claimed"] = True
        with self.assertRaises(ValueError):
            audit.validate(bad)

    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
