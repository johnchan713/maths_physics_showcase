#!/usr/bin/env python3
"""Exact tests for the first-coefficient C3 and leading-C5 boundary."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
import boundary


class FirstCoefficientC3BoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()
        cls.review = cls.record["reviews"]

    def test_derivative_ladder(self):
        self.assertEqual(
            [row["leading_U_order"] for row in self.review["derivative_ladder"]],
            [2, 3, 4, 5, 6, 7],
        )
        self.assertEqual(boundary.derivative_ladder(3)["leading_U_order"], 5)

    def test_reciprocal_and_regular_v_caps(self):
        self.assertEqual(boundary.reciprocal_cap(F(1, 100)), F(25, 23))
        self.assertLess(boundary.regular_v4_bound(F(1, 100), 1)["factor"], 29)

    def test_axis_viscosity_and_weighted_cancellation(self):
        visc = boundary.signed_source_c3_bounds(F(1, 200), 1, 7, 3)
        inviscid = boundary.signed_source_c3_bounds(F(1, 200), 0, 7, 3)
        self.assertGreater(visc["I0_Ceta3"], inviscid["I0_Ceta3"])
        self.assertEqual(visc["I1_Ceta3"], inviscid["I1_Ceta3"])

    def test_first_repair_rows_and_W3(self):
        moments = boundary.first_repair_moments_c3(2, 3, 5, 7, 11, 13)
        self.assertEqual(moments, (6, 12, F(131, 2), 144, F(301, 2)))
        repair = boundary.linear_repair_c3_bound(17, F(1, 100), 19, 23, moments)
        self.assertEqual(repair["W3"], F(5, 23))
        self.assertGreater(repair["coefficient_Ceta3"], 0)

    def test_local_C3_and_cutoff_product(self):
        rows = ((1, 2, 3), (2, 3, 5), (3, 5, 7), (5, 7, 11))
        local = boundary.local_c3_rows(rows, F(1, 4))
        self.assertEqual(local["Ceta2_factor"], 73)
        self.assertEqual(local["Ceta3_factor"], 585)
        with_k3 = boundary.cutoff_product(
            local["radial_Ceta3_bounds"], (1, 729, 2916, 58320)
        )
        without_k3 = boundary.cutoff_product(
            local["radial_Ceta3_bounds"], (1, 729, 2916, 0)
        )
        self.assertEqual(
            with_k3[3][0] - without_k3[3][0],
            F(58320) * local["radial_Ceta3_bounds"][0][0],
        )

    def test_C4_does_not_control_C5(self):
        rows = [boundary.c4_not_c5_countercontrol(n) for n in (1, 8, 64, 512)]
        self.assertTrue(all(row["Ceta4_upper"] <= F(65, 24) for row in rows))
        self.assertEqual(rows[-1]["normalized_fifth_at_zero"], F(64, 15))
        self.assertGreater(
            rows[-1]["normalized_fifth_at_zero"],
            rows[0]["normalized_fifth_at_zero"],
        )

    def test_invalid_inputs(self):
        for value in (-1, True, 1.5):
            with self.assertRaises((ValueError, TypeError)):
                boundary.derivative_ladder(value)
        for args in ((0, 1, 1, 1), (F(1, 200), -1, 1, 1), (F(1, 200), 1, 0, 1)):
            with self.assertRaises(ValueError):
                boundary.signed_source_c3_bounds(*args)
        with self.assertRaises(TypeError):
            boundary.reciprocal_cap(0.01)

    def test_all_exact_checks(self):
        self.assertTrue(all(self.review["checks"].values()))

    def test_scope_not_promoted(self):
        for key in (
            "selected_parameter_whole_support_leading_K05_quantified",
            "finalized_global_first_K13_second_repair_and_H23_closed",
            "raw_deficit_first_band_and_uniform_pulse_family_verified",
            "complete_covariance_axial_pressure_means_residual_or_blowup_verified",
        ):
            self.assertIs(self.record["flags"][key], False)

    def test_dependencies_bound(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)

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
        record["flags"]["selected_parameter_whole_support_leading_K05_quantified"] = 0
        with self.assertRaises(ValueError):
            audit.validate(record)

    def test_protocol_bound(self):
        for key in audit.protocol():
            record = deepcopy(self.record)
            record["protocol"][key] = "erased"
            with self.assertRaises(ValueError):
                audit.validate(record)

    def test_fresh_replay_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory(dir=audit.PROJECT) as directory:
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
