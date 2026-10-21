#!/usr/bin/env python3
"""Tests for the continuation and same-profile matching-root recheck."""

from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import (
    annular_certificate, compact_certificate, continuation_ledger,
    pressure_defect, profile_ledger, review, row_transforms,
)


class CoreContinuationMatchingRecheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_flattened_reference_uses_correct_physical_shear(self):
        row = continuation_ledger()["reference_after_taper"]
        self.assertEqual(row["Dy_log_Fr"], 0)
        self.assertEqual(row["Dy_log_Er"], Q(1, 2))
        self.assertEqual(row["physical_a_r"], 0)
        self.assertFalse(row["first_collar_shortcut_valid_after_taper"])

    def test_five_actual_densities_and_pressure_memory(self):
        data = continuation_ledger()
        self.assertEqual(len(data["moment_densities"]), 5)
        self.assertEqual(data["moment_densities"][-1], "E^2/2")
        self.assertEqual(data["pressure"], "Pi=Pi_axis+Cp")

    def test_final_continuation_slope_matches_shape_entry(self):
        self.assertEqual(continuation_ledger()["shape_join_slope"], Q(1, 10))

    def test_annular_root_is_strictly_contracting(self):
        data = annular_certificate()
        self.assertLess(data["contraction"], 1)
        self.assertEqual(data["radius"], Q(8723, 25) * data["epsilon"])

    def test_annular_root_is_the_old_accepted_root(self):
        data = annular_certificate()
        self.assertLess(data["radius"], data["old_C2_radius"])
        self.assertLess(data["old_C2_radius"], data["old_C1_radius"])

    def test_eta_dependent_annular_quadratic_is_not_frozen(self):
        self.assertEqual(annular_certificate()["Q_eta_normalized_derivatives"], (4, 8, 8, 4, 1))

    def test_compact_root_binary_contraction_margin(self):
        data = compact_certificate()
        self.assertLess(data["contraction_coefficient"], 2**37)
        upper = Q(2**37, data["H_binary_floor"] ** data["H_inverse_power"])
        self.assertLessEqual(upper, Q(1, 2**219))

    def test_compact_transform_needs_positive_lambda(self):
        self.assertTrue(compact_certificate()["lambda_must_be_positive"])
        self.assertIn("lambda", row_transforms()["compact"][1])

    def test_both_maps_restore_the_pressure_row(self):
        rows = row_transforms()
        self.assertEqual(rows["annular"][-1], "Cp")
        self.assertEqual(rows["compact"][-1], "Cp")

    def test_same_profile_reaches_compact_repair(self):
        data = profile_ledger()
        self.assertEqual(len(data["strict_order"]), 6)
        self.assertEqual(data["heat_compensation_starts_at"], "Xrep")
        self.assertIn("[Xminus,Xrep]", data["whole_compact_input"])

    def test_zero_reserved_axial_field_does_not_reset_moment(self):
        data = profile_ledger()
        self.assertEqual(data["reserved_U"], "0")
        self.assertEqual(data["prefix_M"], "mu*eta with mu>0")
        self.assertFalse(data["reserved_M_is_zero"])

    def test_missing_pressure_row_changes_axial_stress(self):
        self.assertEqual(pressure_defect(0, 0, 1), -1)
        self.assertEqual(pressure_defect(1, 1, 0), 2)

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
