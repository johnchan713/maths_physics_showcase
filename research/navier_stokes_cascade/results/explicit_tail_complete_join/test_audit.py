#!/usr/bin/env python3
"""Tests for the explicit-tail complete canonical join."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from join_checks import B0, KAPPA, combined_tail, cycle_orders, force_terms, review


class ExplicitTailCompleteJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_combined_tail_is_maximum_not_object_identification(self):
        self.assertEqual(combined_tail(11, 31, 17), 31)
        self.assertEqual(combined_tail(10**100, 2, 3), 10**100)
        self.assertFalse(audit.FLAGS["inherited_and_canonical_primary_fields_or_roots_identified"])

    def test_combined_tail_input_contract(self):
        with self.assertRaises(TypeError):
            combined_tail()
        with self.assertRaises(TypeError):
            combined_tail(3.0, 4)
        with self.assertRaises(ValueError):
            combined_tail(1, 4)

    def test_cycle_seed_and_all_groups_close(self):
        row = cycle_orders(B0)
        self.assertTrue(all(margin >= 0 for margin in row["margins"].values()))
        self.assertEqual(row["tie_count"], 7)
        self.assertEqual(row["orders"]["new_square"] - row["required"], F(1) - 2 * KAPPA)

    def test_cycle_affine_half_line(self):
        for j in (0, 1, 2, 10, 100, 10**5):
            B = B0 + j * (F(1, 2) - KAPPA)
            row = cycle_orders(B)
            self.assertTrue(all(margin >= 0 for margin in row["margins"].values()))
            self.assertEqual(row["tie_count"], 7)

    def test_cycle_order_rejects_float_and_below_seed(self):
        with self.assertRaises(TypeError):
            cycle_orders(float(B0))
        with self.assertRaises(ValueError):
            cycle_orders(B0 - F(1, 10**9))

    def test_localization_ratio_and_force_ledger(self):
        exact = review()
        self.assertEqual(exact["localization"]["q_over_qsharp"], F(3, 128))
        self.assertEqual(exact["localization"]["whole_chart_over_qsharp"], F(3, 8))
        self.assertLess(exact["localization"]["whole_chart_over_qsharp"], 1)
        self.assertEqual(len(force_terms()), 12)
        self.assertIn("p_gradient_c", force_terms())
        self.assertIn("laplacian_d", force_terms())

    def test_dependency_routes_and_source_bindings(self):
        self.assertTrue(all(
            self.record["reviews"]["dependency_routing"]["checks"].values()
        ))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 11)

    def test_scope_remains_conditional(self):
        for key in (
            "inherited_and_canonical_primary_fields_or_roots_identified",
            "historical_first_band_realization_force_and_schedules_reidentified",
            "fixed_stage_estimate_promoted_to_uncontrolled_stage_diagonal",
            "retained_core_background_defining_map_and_inverse_premises_all_discharged",
            "unconditional_pulse_complete_residual_admissible_blowup_or_blowup_proved",
            "formal_certificate",
        ):
            self.assertIs(self.record["flags"][key], False)

    def test_false_scope_promotions_rejected(self):
        for key, value in audit.FLAGS.items():
            if not value:
                changed = deepcopy(self.record)
                changed["flags"][key] = True
                with self.assertRaises(ValueError):
                    audit.validate(changed)

    def test_pressure_or_force_ledger_tamper_rejected(self):
        changed = deepcopy(self.record)
        changed["reviews"]["exact_join"]["force_terms"].remove("p_gradient_c")
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_remaining_foundation_obstruction_cannot_be_erased(self):
        changed = deepcopy(self.record)
        changed["remaining_obligations"].pop(0)
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
