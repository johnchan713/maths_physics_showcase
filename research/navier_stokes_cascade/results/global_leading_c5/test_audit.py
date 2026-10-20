#!/usr/bin/env python3
"""Tests for the selected whole-support leading C_eta^5 checkpoint."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
import leading5


class GlobalLeadingC5Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()
        cls.review = cls.record["reviews"]

    def test_exact_fifth_power_chain_rule(self):
        self.assertEqual(
            leading5.fifth_power_terms(),
            ((5, 0, 5, 32), (3, 1, 4, 160), (1, 2, 3, 120)),
        )
        self.assertEqual(self.review["later_reference"]["raw_derivative_caps"][-1], 3936)
        self.assertLess(F(self.review["later_reference"]["interpolation_Ceta5"]), 72)

    def test_same_strength_root_five_jets(self):
        eps = F(1, 2**20)
        jets = leading5.strength_jets(eps)
        self.assertEqual(len(jets), 6)
        self.assertTrue(all(jets[n] <= 12 * eps for n in range(1, 6)))

    def test_invalid_strength_input(self):
        for value in (-1, F(1, 2**19), True, 0.5):
            with self.assertRaises((ValueError, TypeError)):
                leading5.strength_jets(value)

    def test_complete_later_reference_cover(self):
        later = self.review["later_reference"]
        self.assertEqual(later["later_reference_U_Ceta5"], 132)
        self.assertEqual(later["complete_energy_Ceta5"], "8192 T")
        self.assertIn("SAME", later["angular_root_Ceta5"])

    def test_whole_support_join(self):
        joined = self.review["whole_support"]
        self.assertEqual(
            joined["whole_support_K05"],
            "sup_[0,Xv] |U0|Ceta5 < H^262",
        )
        self.assertIn("H^261", joined["FINAL_compact_U_Ceta5"])

    def test_all_exact_checks(self):
        self.assertTrue(all(self.review["checks"].values()))

    def test_scope_not_promoted(self):
        for key in (
            "finalized_global_first_K13_closed",
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
                record = deepcopy(self.record)
                record["flags"][key] = True
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
