#!/usr/bin/env python3
"""Tests for the independent core growth-source recheck."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import exit_certificate, review, scalar_certificate, source_certificate


class CoreGrowthSourceRecheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_scalar_lower_and_log_slope_are_exact(self):
        data = scalar_certificate()
        self.assertGreater(data["f0_lower_cubic"], F(53, 200))
        self.assertLess(data["log_slope_upper"], 4)
        self.assertEqual(data["weighted_eta_upper"], 410)

    def test_endpoint_alternating_polynomial(self):
        data = scalar_certificate()
        self.assertEqual(data["endpoint_polynomial"], -F(75535511, 400000000))
        self.assertGreater(data["endpoint_comparison_shear_lower"], F(2376, 1000))

    def test_complete_source_group_ledger(self):
        data = source_certificate()
        self.assertEqual(len(data["terms"]), 10)
        self.assertEqual(len({row[0] for row in data["terms"]}), 10)
        self.assertEqual(data["chi_loss"], 520)
        self.assertEqual(data["fixed_point_error_multiplier"], 2160)

    def test_eta_weight_and_reciprocal_error(self):
        data = scalar_certificate()
        self.assertEqual(data["weighted_reciprocal_error"], 1896)
        self.assertLess(data["weighted_reciprocal_error"], 2000)
        self.assertTrue(self.record["flags"]["weighted_eta_reciprocal_bound_keeps_Hstar_factor"])

    def test_young_and_universal_source_remainders(self):
        data = source_certificate()
        self.assertLess(data["young_remainder_coefficient"], 9)
        self.assertLess(data["universal_remainder_prefactor"], F(3, 1000))

    def test_forward_integral_and_native_a(self):
        data = exit_certificate()
        self.assertEqual(data["p1_over_X_lower"], F(5, 32))
        self.assertGreater(data["p1_over_X_lower"], F(1, 10))
        self.assertEqual(data["native_a_lower_coefficient_over_Lambda"], F(2, 5))

    def test_axial_physical_scaling_denominator(self):
        data = exit_certificate()
        self.assertEqual(data["axial_ratio_denominator"], 80)
        self.assertEqual(data["native_v_lower"], F(23, 10))

    def test_every_exact_gate_passes(self):
        self.assertTrue(all(review()["checks"].values()))
        self.assertTrue(all(self.record["gates"].values()))

    def test_false_scope_promotion_is_rejected(self):
        for key, value in audit.FLAGS.items():
            if not value:
                changed = deepcopy(self.record)
                changed["flags"][key] = True
                with self.assertRaises(ValueError):
                    audit.validate(changed)

    def test_remaining_frontier_cannot_be_erased(self):
        changed = deepcopy(self.record)
        changed["remaining_obligations"].pop()
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_frozen_record_and_fresh_replay(self):
        audit.validate(json.loads((audit.HERE / "evidence.json").read_text()))
        with tempfile.TemporaryDirectory(dir=audit.PROJECT) as directory:
            output = Path(directory) / "fresh.json"
            audit.main([
                "--verify-record", str(audit.HERE / "evidence.json"),
                "--output", str(output),
            ])
            audit.validate(json.loads(output.read_text()))
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(output)])


if __name__ == "__main__":
    unittest.main()
