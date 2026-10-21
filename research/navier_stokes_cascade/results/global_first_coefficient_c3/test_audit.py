#!/usr/bin/env python3
"""Tests for the finalized whole-support first U coefficient C_eta^3 bound."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
import global_first


class GlobalFirstCoefficientC3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()
        cls.review = cls.record["reviews"]

    def test_local_value_charges_only_angular_cauchy_and_cutoff_value(self):
        local = global_first.local_value_c3(F(1, 4), 7, 11)
        self.assertEqual(local["Ceta3_factor"], 585)
        self.assertEqual(local["cutoff_value_cost"], 1)
        self.assertEqual(local["tentative_value_Ceta3"], 585 * 18)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            global_first.local_value_c3(0, 1, 1)
        with self.assertRaises(TypeError):
            global_first.local_value_c3(0.25, 1, 1)
        with self.assertRaises(ValueError):
            global_first.inverse_800(1, 1, 1, 1, (1, 2, 3, 4))

    def test_full_source_retains_axis_viscosity_only_in_unweighted_total(self):
        full = global_first.signed_source_c3(F(1, 200), 1, 13, 17)
        inviscid = global_first.signed_source_c3(F(1, 200), 0, 13, 17)
        self.assertGreater(full["Q03"], inviscid["Q03"])
        self.assertEqual(full["Q13"], inviscid["Q13"])

    def test_actual_bump_value_and_linear_repair(self):
        out = self.review["exact_fixture"]["assembled"]
        local = self.review["exact_fixture"]["local"]
        self.assertEqual(out["actual_U_bump_value"], "320")
        self.assertEqual(
            F(out["K13"]),
            F(local["tentative_value_Ceta3"])
            + 320 * F(out["repair"]["coefficient_Ceta3"]),
        )
        self.assertEqual(F(out["repair"]["W3"]), F(5, 41))

    def test_selected_formula_keeps_small_denominators_visible(self):
        selected = self.review["selected_formula"]
        self.assertEqual(selected["leading_inputs"]["K05"], "H^262")
        self.assertIn("lambda", selected["repair"])
        self.assertEqual(selected["final"], "K1,3 := J13 + 320 C13")

    def test_all_exact_checks(self):
        self.assertTrue(all(self.review["checks"].values()))

    def test_scope_not_promoted(self):
        for key in (
            "second_signed_moments_repair_and_H23_closed",
            "uniform_actual_physical_pulse_family_verified",
            "complete_covariance_means_residual_or_blowup_verified",
        ):
            self.assertIs(self.record["flags"][key], False)

    def test_dependencies_bound(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 7)

    def test_frozen_record(self):
        audit.validate(json.loads((audit.HERE / "evidence.json").read_text()))

    def test_false_scope_promotions_rejected(self):
        for key, value in audit.FLAGS.items():
            if not value:
                changed = deepcopy(self.record)
                changed["flags"][key] = True
                with self.assertRaises(ValueError):
                    audit.validate(changed)

    def test_protocol_bound(self):
        for key in audit.protocol():
            changed = deepcopy(self.record)
            changed["protocol"][key] = "erased"
            with self.assertRaises(ValueError):
                audit.validate(changed)

    def test_fresh_replay_and_overwrite_refusal(self):
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
