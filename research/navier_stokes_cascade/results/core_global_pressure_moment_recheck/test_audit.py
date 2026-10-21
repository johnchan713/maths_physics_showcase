#!/usr/bin/env python3
"""Tests for the same-profile global pressure and moment recheck."""

from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import (
    annular_inverse, annular_rows, compact_inverse, compact_rows,
    edit_ledger, endpoint_to_global_ledger, exterior_defects,
    heat_ledger, pressure_fixture, review,
)


class CoreGlobalPressureMomentRecheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_global_heat_object_is_renormalized(self):
        data = heat_ledger()
        self.assertFalse(data["raw_I_infinity_exists"])
        self.assertIn("I(X)-", data["R_limit"])
        self.assertEqual(data["axis_Hpow_integrable_condition"], "0<h<1")

    def test_h_zero_limit_is_not_interchanged(self):
        self.assertFalse(heat_ledger()["h_zero_interchange_allowed"])

    def test_endpoint_I_and_global_R_are_distinct_rows(self):
        data = endpoint_to_global_ledger()
        self.assertEqual(data["endpoint_row_order"], "(M,I,J,S,Cp)")
        self.assertEqual(data["global_defect_order"], "(M,R,J,S,Cp)")
        self.assertTrue(data["common_future_tail"])
        self.assertTrue(data["common_Hpow"])

    def test_raw_I_infinity_is_not_zero(self):
        data = endpoint_to_global_ledger()
        self.assertNotEqual(data["invalid_conclusion"], data["conclusion"])

    def test_annular_endpoint_rows_invert_exactly(self):
        values = (Q(2), Q(3), Q(5), Q(7), Q(11))
        self.assertEqual(annular_inverse(annular_rows(values, Q(2, 5)), Q(2, 5)), values)

    def test_compact_endpoint_rows_invert_exactly(self):
        values = (Q(2), Q(3), Q(5), Q(7), Q(11))
        lam = Q(1, 10**20)
        self.assertEqual(compact_inverse(compact_rows(values, lam), lam), values)
        with self.assertRaises(ValueError):
            compact_rows(values, 0)
        with self.assertRaises(ValueError):
            compact_inverse(values, 0)

    def test_five_global_defects_are_functional_equalities(self):
        data = edit_ledger()
        self.assertEqual(data["global_zero_defects"], ("M", "R", "J", "S", "Pi_axis+Cp"))
        self.assertTrue(data["functional_in_eta"])

    def test_same_pressure_gauge_gives_backward_tail(self):
        data = pressure_fixture()
        self.assertEqual(data["pi_axis"] + data["final_total"], 0)
        self.assertEqual(data["forward"], data["backward"])

    def test_pressure_memory_is_not_reset_between_supports(self):
        self.assertTrue(any(pressure_fixture()["intermediate_heat"]))

    def test_omitting_pressure_repair_breaks_gauge(self):
        self.assertNotEqual(pressure_fixture()["missing_pressure_repair"], 0)

    def test_zero_defects_remove_all_exterior_memory(self):
        zero = [(Q(0), Q(0))] * 5
        self.assertEqual(exterior_defects(Q(2, 5), Q(1, 1000), 100, Q(3, 2), zero), (0, 0, 0))

    def test_pressure_slope_defect_grows_linearly(self):
        defects = [(Q(0), Q(0))] * 5
        defects[4] = (Q(0), Q(1))
        z100 = exterior_defects(0, Q(1, 1000), 100, 1, defects)[2]
        z200 = exterior_defects(0, Q(1, 1000), 200, 1, defects)[2]
        self.assertNotEqual(z100, 0)
        self.assertEqual(z200, 2 * z100)

    def test_zero_U_does_not_erase_M_jet(self):
        defects = [(Q(0), Q(0))] * 5
        defects[0] = (Q(0), Q(1))
        self.assertNotEqual(exterior_defects(0, Q(1, 1000), 100, 1, defects)[0], 0)

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
