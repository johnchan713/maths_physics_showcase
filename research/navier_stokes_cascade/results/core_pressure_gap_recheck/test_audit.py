#!/usr/bin/env python3
"""Tests for the independent shared-core pressure and gap recheck."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import (
    A_MAX, J_MAX, atom_pressure_certificate, complex_pressure_certificate,
    direct_velocity_polynomial, parameter_order, review,
    separation_certificate, stated_velocity_polynomial, velocity_envelope,
)


class CorePressureGapRecheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_velocity_polynomial_is_independently_expanded(self):
        for A in (F(1, 2), A_MAX, F(503, 1000)):
            for j in (F(0), F(1, 10**18), J_MAX):
                self.assertEqual(direct_velocity_polynomial(A, j),
                                 stated_velocity_polynomial(A, j))

    def test_sharper_velocity_envelope(self):
        bounds = velocity_envelope()
        self.assertLess(bounds["eta_total"], 21)
        self.assertLess(bounds["constant_over_j"], 5)

    def test_atom_pressure_minimum_is_endpoint(self):
        atom = atom_pressure_certificate()
        self.assertLess(atom["derivative_numerator_max_without_factor_10"], 0)
        self.assertEqual(atom["endpoint_min"], F(5, 4))

    def test_sharper_separation_and_chi_floor(self):
        gap = separation_certificate()
        self.assertEqual(gap["exact_ratio"], F(189, 7475))
        self.assertLess(gap["exact_ratio"], F(1, 39))
        self.assertGreater(gap["H_over_j_exact_lower_control"], F(22, 25))
        self.assertEqual(gap["chi_floor"], F(7744, 7745))
        self.assertGreater(gap["chi_floor"], F(99, 100))

    def test_complex_pressure_and_parameter_order(self):
        pressure = complex_pressure_certificate()
        self.assertLess(pressure["mixture_upper_over_P2"], 10)
        order = parameter_order()
        self.assertTrue(order["all_edges_forward"])
        self.assertEqual(order["radial_scale_power_in_dX_over_X"], 0)

    def test_every_exact_check_passes(self):
        self.assertTrue(all(review()["checks"].values()))
        self.assertTrue(all(self.record["gates"].values()))

    def test_scope_promotion_is_rejected(self):
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
