#!/usr/bin/env python3
"""Exact tests for the explicit coefficient-cutoff selection."""

from copy import deepcopy
import json
import tempfile
import unittest
from pathlib import Path

import audit
import cutoff


class ExplicitCoefficientCutoffTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()
        cls.selected = cls.record["reviews"]["selected_coefficient_cutoff"]

    def test_exact_bump_majorants(self):
        self.assertEqual(
            [row["integer_majorant"] for row in self.selected["bump_rows"]],
            [
                1,
                4,
                80,
                3600,
                292320,
                37704960,
                7139819520,
                1867003084800,
            ],
        )

    def test_normalized_cutoff_jets(self):
        self.assertEqual(
            self.selected["normalized_derivative_strict_upper"][:4],
            [1, 729, 2916, 58320],
        )

    def test_stirling_chain_rule(self):
        self.assertEqual(
            self.selected["signed_stirling_rows"][3], [0, 2, -3, 1]
        )
        self.assertEqual(
            self.selected["physical_derivative_upper_polynomials"][3]["terms"],
            [
                {"delta_power": -1, "coefficient": 1458},
                {"delta_power": -2, "coefficient": 8748},
                {"delta_power": -3, "coefficient": 58320},
            ],
        )

    def test_all_exact_checks(self):
        self.assertTrue(all(self.selected["checks"].values()))

    def test_invalid_orders(self):
        for order in (True, 2, 9, 3.0):
            with self.assertRaises(ValueError):
                cutoff.review(order)

    def test_cutoffs_are_distinguished(self):
        self.assertIs(
            self.record["flags"]["coefficient_and_variance_cutoffs_are_not_conflated"],
            True,
        )

    def test_precise_remaining_gate(self):
        gate = self.record["reviews"]["first_unclosed_mixed_gate"]
        self.assertEqual((gate["M"], gate["N"]), (729, 2916))
        self.assertIn("d0", gate["budget"])
        self.assertIn("H>0", gate["required"])

    def test_no_scope_promotion(self):
        for key in (
            "actual_first_scale_raw_deficit_and_mixed_gate_verified",
            "global_first_U1_Ceta3_and_finalized_second_H23_verified",
            "uniform_complete_pulse_family_covariance_means_residual_or_blowup_verified",
        ):
            self.assertIs(self.record["flags"][key], False)

    def test_nine_dependency_records_bound(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 9)

    def test_frozen_record(self):
        audit.validate(json.loads((audit.HERE / "evidence.json").read_text()))

    def test_false_scope_promotions_rejected(self):
        for key, value in audit.FLAGS.items():
            if not value:
                record = deepcopy(self.record)
                record["flags"][key] = True
                with self.assertRaises(ValueError):
                    audit.validate(record)

    def test_zero_cannot_replace_false(self):
        record = deepcopy(self.record)
        record["flags"][
            "actual_first_scale_raw_deficit_and_mixed_gate_verified"
        ] = 0
        with self.assertRaises(ValueError):
            audit.validate(record)

    def test_protocol_is_bound(self):
        for key in audit.protocol():
            record = deepcopy(self.record)
            record["protocol"][key] = "erased"
            with self.assertRaises(ValueError):
                audit.validate(record)

    def test_fresh_replay_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "fresh.json"
            audit.main(
                [
                    "--verify-record",
                    str(audit.HERE / "evidence.json"),
                    "--output",
                    str(output),
                ]
            )
            audit.validate(json.loads(output.read_text()))
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(output)])


if __name__ == "__main__":
    unittest.main()
