#!/usr/bin/env python3
"""Tests for the full-support second repair and explicit global H2,3."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
import global_second


class GlobalSecondCoefficientH23Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()
        cls.review = cls.record["reviews"]

    def test_regular_average_requires_extra_leading_jet(self):
        self.assertGreater(
            global_second.regular_v_bound(0, F(1, 200), 3, 1),
            global_second.regular_v_bound(0, F(1, 200), 2, 1),
        )

    def test_both_shifted_viscosity_moments_are_nonzero(self):
        full = global_second.shifted_second_source(F(1, 200), 1, 7, 11, 13)
        inviscid = global_second.shifted_second_source(F(1, 200), 0, 7, 11, 13)
        self.assertGreater(full["Q10"], inviscid["Q10"])
        self.assertGreater(full["Q11"], inviscid["Q11"])
        self.assertTrue(all(value > 0 for value in full["shifted_axial_viscosity"]))

    def test_first_F_keeps_inverse_radius_cost(self):
        fields = global_second.first_field_values(2, 3, F(1, 2))
        self.assertEqual(fields["U1"], 962)
        self.assertEqual(fields["F1"], 1922)

    def test_cutoff_third_jet_keeps_log_jacobian_terms(self):
        jets = global_second.coefficient_cutoff_jets(F(1, 8), 2)
        top_only = F(58320 * 8**3, 2**3)
        self.assertGreater(jets[3], top_only)

    def test_product_keeps_k3_J20(self):
        cut = global_second.product_jets((1, 2, 3, 5), (7, 11, 13, 17))
        omit = global_second.product_jets((1, 2, 3, 0), (7, 11, 13, 17))
        self.assertEqual(cut[3] - omit[3], 35)

    def test_second_collar_is_evaluated_from_finalized_first_rows(self):
        first = global_second.finalized_first_jets(
            (2, 3, 5, 7), (1, 2, 3, 5), 11, 13, F(1, 2)
        )
        t_collar = max(F(1), 17, first["mixed_Ceta3"])
        second = global_second.second_inner_jets(
            (2, 3, 5, 7), ((11, 13), (17, 19), (23, 29), (31, 37)),
            t_collar, 41, 43,
        )
        self.assertEqual(t_collar, first["mixed_Ceta3"])
        self.assertEqual(second[0], 2 * 41 + 11 * 43 + 13)

    def test_lower_first_squares_and_products_remain(self):
        rows = global_second.second_discrepancies(2, 3, 5, 7, 11, 13, 17, 19, 23)
        self.assertGreater(rows[2], 2 * 7 * 2 * 5 + F(19, 2))
        self.assertGreater(rows[4], 2 * 11 * 2 * 5 + 7 * 4 * 5 + F(23, 2))

    def test_linear_inverse_uses_W1(self):
        out = global_second.inverse_800(2, F(1, 1000), 3, 5, (7, 11, 13, 17, 19))
        self.assertEqual(out["W1"], F(4, 5))
        self.assertEqual(out["coefficient_Ceta1"], 800 * max(out["rows"]))

    def test_actual_third_bump_and_E_to_F_cost(self):
        jets = global_second.repair_bump_jets(2, F(1, 2))
        self.assertEqual(jets["logarithmic"][3], 18870643200000)
        self.assertGreater(jets["E_to_F"][3], jets["U"][3])

    def test_final_H23_converts_F_to_E(self):
        converted = global_second.f_to_e_jets((2, 3, 5, 7), F(1, 2), 3)
        self.assertEqual(converted["E"][0], 6)
        self.assertEqual(converted["radius_jets"], (3, 2, 8, 96))

    def test_all_exact_checks(self):
        self.assertTrue(all(self.review["checks"].values()))

    def test_scope_not_promoted(self):
        for key in (
            "historical_first_physical_band_raw_deficit_gate_verified",
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

    def test_invalid_inputs(self):
        with self.assertRaises(TypeError):
            global_second.regular_v_bound(0, 0.005, 2, 1)
        with self.assertRaises(ValueError):
            global_second.inverse_800(1, 1, 1, 1, (1, 2, 3, 4))
        with self.assertRaises(ValueError):
            global_second.coefficient_cutoff_jets(0, 1)

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
