#!/usr/bin/env python3
"""Tests for the same-profile modulation, repair, and cone-join recheck."""

from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import (
    compact_inverse, compact_transform, fast_chain_control,
    hierarchy_certificate, moving_phase_control, pressure_defect,
    repair_certificate, review, signed_variance_ledger, support_ledger,
)


class CoreModulationRepairJoinRecheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_signed_root_is_regular_at_zero(self):
        data = signed_variance_ledger()
        self.assertEqual(data["R2_at_zero"], Q(1, 2))
        self.assertEqual(data["signed_F_mu_squared_over_d0_squared_at_zero"], Q(1, 2))
        self.assertFalse(data["division_by_V_mu_at_zero_allowed"])

    def test_phase_inverse_must_move(self):
        data = moving_phase_control()
        self.assertEqual(data["inverse_identity"], 0)
        self.assertNotEqual(data["phase_term"], 0)
        self.assertNotEqual(data["pulled_t_r"], data["t_r_fixed_theta"])

    def test_fast_chain_term_is_retained(self):
        data = fast_chain_control(7)
        self.assertEqual(data["full_radial_derivative"], Q(15, 7))
        self.assertEqual(data["slow_remainder"], Q(1, 7))
        self.assertNotEqual(data["full_radial_derivative"], data["slow_remainder"])

    def test_moments_gain_one_radial_order(self):
        data = fast_chain_control()
        self.assertEqual(data["field_rth_radial_power"], "N^(r-1)")
        self.assertEqual(data["moment_rth_radial_power"], "N^max(r-2,-1)")

    def test_sharper_same_root_contraction(self):
        data = repair_certificate()
        self.assertEqual(data["root_coefficient"], 1600)
        self.assertLess(data["contraction_coefficient"], 2**37)
        self.assertLessEqual(Q(2**37, (2**16)**16), Q(1, 2**219))

    def test_compact_rows_invert_for_positive_lambda(self):
        values = (Q(2), Q(3), Q(5), Q(7), Q(11))
        rows = compact_transform(values, Q(1, 100))
        self.assertEqual(compact_inverse(rows, Q(1, 100)), values)

    def test_compact_rows_reject_zero_lambda(self):
        with self.assertRaises(ValueError):
            compact_transform((2, 3, 5, 7, 11), 0)
        with self.assertRaises(ValueError):
            compact_inverse((2, 3, 5, 7, 11), 0)

    def test_pressure_row_is_not_optional(self):
        self.assertEqual(repair_certificate()["rows"][-1], "Cp")
        self.assertEqual(pressure_defect(0, 0, 1), -1)

    def test_selected_hierarchy_clears_state_and_radius_thresholds(self):
        data = hierarchy_certificate()
        self.assertLess(data["state_polynomial_left"], data["state_polynomial_right"])
        self.assertLess(data["radius_polynomial_left"], data["radius_polynomial_right"])

    def test_cone_constant_is_absorbed_at_floor(self):
        data = hierarchy_certificate()
        self.assertGreater(data["cone_left_exponent"], data["cone_constant_exponent"])
        self.assertEqual(data["cone_constant_exponent"], 29)

    def test_supports_join_without_reselection(self):
        data = support_ledger()
        self.assertEqual(data["repair_endpoint"], "Xrep=rho_p*exp(5)=Xstar")
        self.assertFalse(data["annular_and_compact_roots_reselected"])
        self.assertFalse(data["frequency_reselected"])
        self.assertEqual(data["inner_width"], "C^-20")

    def test_all_five_endpoint_defects_are_functional_equalities(self):
        data = support_ledger()
        self.assertTrue(data["all_five_endpoint_defects_zero_as_eta_functions"])
        self.assertTrue(data["axis_pressure_datum_unchanged"])

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
