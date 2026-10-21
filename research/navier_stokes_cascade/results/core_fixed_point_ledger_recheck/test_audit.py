#!/usr/bin/env python3
"""Tests for the independent all-index core fixed-point ledger."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import (
    A0, cauchy_multiplier, convolution_sum, group_sums, inverse_word_bound,
    ledger, lipschitz_coefficients, majorant_coefficients, mixed_allocation,
    radial_shift_ratio, review, universal_contraction,
)


class CoreFixedPointLedgerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_all_source_terms_are_present_once(self):
        rows = ledger()
        self.assertEqual(len(rows), 21)
        self.assertEqual(len({row["name"] for row in rows}), 21)

    def test_exact_source_group_sums(self):
        sums = group_sums()
        self.assertEqual(sums["R1_linear"], F(19237, 50))
        self.assertEqual(sums["R1_quadratic"], F(952, 25))
        self.assertEqual(sums["R1_zeta"], 10)
        self.assertEqual(sums["R2_linear"], 442)
        self.assertEqual(sums["R2_quadratic"], 32)
        self.assertEqual(sums["R2_pressure"], 34)

    def test_pressure_algebra_and_lipschitz_factors(self):
        major = majorant_coefficients()
        lip = lipschitz_coefficients()
        self.assertEqual(major["pressure_primitive_A0_cubed"], 80 * A0**3)
        self.assertEqual(major["R2_pressure_A0_fourth"], 2720 * A0**4)
        self.assertEqual(lip["R2_pressure_per_S"], 5440 * A0**4)

    def test_omitting_each_pressure_term_is_detected(self):
        rows = ledger()
        pressure = [row for row in rows if row["group"] == "R2_pressure"]
        self.assertEqual({row["constant"] for row in pressure}, {F(16), F(10), F(8)})
        for row in pressure:
            reduced = tuple(item for item in rows if item["name"] != row["name"])
            self.assertNotEqual(group_sums(reduced)["R2_pressure"], 34)

    def test_weight_shifts_and_integrated_derivative_allocation(self):
        for n in range(16):
            for k in range(16):
                self.assertEqual(radial_shift_ratio(n, k),
                                 F(20 * (n + 2)**2, (n + 1) * (n + k + 1)))
                self.assertLessEqual(radial_shift_ratio(n, k), 80)
        for i in (0, 1, 9, 10**6):
            for j in (0, 1, 13, 10**6):
                for nu in (1, 2):
                    self.assertLessEqual(mixed_allocation(i, j, nu), 1)

    def test_product_and_cauchy_controls(self):
        for n in range(65):
            self.assertLess(convolution_sum(n), 16)
        for k in range(65):
            self.assertLessEqual(cauchy_multiplier(k), 1)

    def test_inverse_word_factorials(self):
        for k in range(20):
            self.assertEqual(inverse_word_bound(k + 1) / inverse_word_bound(k),
                             F(30720, (k + 1) * (k + 2)))

    def test_selected_formula_has_strict_universal_ceiling(self):
        result = universal_contraction()
        self.assertEqual(result["epsilon_squared_over_10_power_18"], F(1, 10**50))
        self.assertEqual(result["H_denominator_power"], 127)

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
