from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from slack import (
    chosen_log_scale, full_exponent_model, review, scale_margin, slack_gain,
    slack_model, stage_gain, young_log_majorant,
)


class CornerTailSlackReconciliationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_positive_gain_is_split_exactly_in_half(self):
        for stage in (1, 4, 17):
            self.assertGreater(stage_gain(stage), 0)
            self.assertEqual(2 * slack_gain(stage), stage_gain(stage))

    def test_full_exponent_no_log_normalization_is_unbounded(self):
        for power in (1, 3, 20):
            self.assertGreater(full_exponent_model(10**6, power),
                               full_exponent_model(10**3, power))

    def test_half_gain_has_a_finite_analytic_majorant(self):
        for stage, power in ((1, 0), (1, 7), (17, 2000)):
            self.assertEqual(young_log_majorant(stage, power),
                             Q(power * power, 1) / stage_gain(stage))

    def test_whole_interval_scale_rule_absorbs_extreme_fixed_data(self):
        previous = Q(0)
        for stage, constant, power in ((1, 5, 7), (4, 1000, 300), (17, 10**8, 2000)):
            chosen = chosen_log_scale(stage, constant, power, previous)
            self.assertLess(scale_margin(stage, constant, power, chosen), 0)
            self.assertGreater(chosen, previous + 1)
            previous = chosen

    def test_slack_model_decays_after_its_fixed_peak(self):
        self.assertLess(slack_model(10**5, 1, 3), slack_model(10**4, 1, 3))

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError): stage_gain(0)
        with self.assertRaises(ValueError): full_exponent_model(0, 1)
        with self.assertRaises(ValueError): young_log_majorant(1, -1)
        with self.assertRaises(ValueError): chosen_log_scale(1, -1, 0)

    def test_controls_dependencies_and_scientific_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 7)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))
        self.assertFalse(result["scope"]["actual_primitive_closed_coefficients_proved"])
        self.assertFalse(result["scope"]["full_exponent_no_log_tail_supremum_required"])

    def test_frozen_record_reproduces_and_scope_tampering_fails(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["actual_primitive_closed_coefficients_proved"] = True
        with self.assertRaises(ValueError): audit.validate(bad)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit): audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
