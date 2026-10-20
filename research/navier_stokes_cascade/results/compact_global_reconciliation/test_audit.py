#!/usr/bin/env python3
"""Tests for the compact/global reconciliation checkpoint."""
import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

import audit
import join


class CompactGlobalReconciliationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()
        cls.review = cls.record["reviews"]

    def test_all_gates(self):
        self.assertTrue(all(self.record["gates"].values()))

    def test_fifth_gate_implies_F5(self):
        self.assertTrue(self.review["checks"]["fifth_gate_implies_F5"])
        self.assertEqual(self.review["bounds"]["compact_ordinary_C4_upper"], "H^925")

    def test_complete_compact_target(self):
        self.assertTrue(self.review["checks"]["complete_compact_C4_H925_below_H1024"])

    def test_cone_margin(self):
        self.assertTrue(self.review["checks"]["cone_perturbation_below_half_gap"])
        self.assertEqual(self.review["bounds"]["cone_gap_lower"], "G^-4/2")

    def test_both_transforms_invert(self):
        self.assertTrue(self.review["checks"]["annular_transform_exactly_invertible"])
        self.assertTrue(self.review["checks"]["compact_transform_exactly_invertible_for_positive_lambda"])

    def test_pressure_negative_control(self):
        self.assertEqual(self.review["identities"]["pressure_defect_control_at_eta0_X13"], "-13")
        self.assertTrue(self.review["checks"]["pressure_defect_negative_control_nonzero"])

    def test_dependency_records_are_bound(self):
        self.assertEqual(set(self.record["inherited_provenance"]["statuses"]), set(audit.DEPENDENCIES))
        self.assertEqual(len(self.record["inherited_provenance"]["frozen_records"]), len(audit.DEPENDENCIES))

    def test_scope_not_promoted(self):
        for key in ("independent_core_native_regional_or_root_foundation_review",
                    "global_FINAL_all_order_leading_profile_or_all_physical_scale_join",
                    "uniform_all_family_pulses_complete_covariance_axial_pressure_means_residual_or_blowup"):
            self.assertIs(self.record["flags"][key], False)

    def test_false_scope_promotion_rejected(self):
        for key, value in audit.FLAGS.items():
            if not value:
                changed = deepcopy(self.record)
                changed["flags"][key] = True
                with self.assertRaises(ValueError):
                    audit.validate(changed)

    def test_pressure_row_cannot_be_removed(self):
        changed = deepcopy(self.record)
        changed["flags"]["all_five_global_moment_rows_and_complete_pressure_memory_retained"] = False
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_protocol_is_bound(self):
        for key in audit.protocol():
            changed = deepcopy(self.record)
            changed["protocol"][key] = "erased"
            with self.assertRaises(ValueError):
                audit.validate(changed)

    def test_remaining_obligations_are_bound(self):
        changed = deepcopy(self.record)
        changed["remaining_obligations"].pop()
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_frozen_evidence(self):
        audit.validate(json.loads((audit.HERE / "evidence.json").read_text()))

    def test_fresh_replay_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "fresh.json"
            audit.main(["--verify-record", str(audit.HERE / "evidence.json"), "--output", str(out)])
            audit.validate(json.loads(out.read_text()))
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(out)])

    def test_review_is_deterministic(self):
        self.assertEqual(join.review(), join.review())


if __name__ == "__main__":
    unittest.main()
