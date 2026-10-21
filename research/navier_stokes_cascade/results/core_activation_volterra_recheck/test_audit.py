#!/usr/bin/env python3
"""Tests for the independent flat-activation and five-moment recheck."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import (
    cone_ledger, envelope_ledger, map_ledger, moment_ledger, review,
    scale_certificate, volterra_certificate,
)


class CoreActivationVolterraRecheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_positive_Volterra_bound_and_endpoint(self):
        data = volterra_certificate()
        self.assertEqual(data["operator_norm"], F(3, 2))
        self.assertEqual(data["endpoint_value_multiplier"], F(1, 2))
        self.assertEqual(data["V0_second_jet_b0_multiplier"], -F(3, 2))

    def test_all_five_moments_and_exact_C_weights(self):
        rows = moment_ledger()
        self.assertEqual(tuple(row["name"] for row in rows), ("M", "I", "J", "S", "Cp"))
        self.assertEqual(tuple(row["C_power"] for row in rows), (0, -1, -1, 0, -2))

    def test_product_and_pressure_moment_integrands(self):
        rows = {row["name"]: row for row in moment_ledger()}
        self.assertEqual(len(rows["J"]["formula"]), 3)
        self.assertEqual(rows["Cp"]["formula"], ("bE2/2",))

    def test_complete_Qs_reciprocal_difference(self):
        rows = map_ledger()
        self.assertEqual(len(rows["Delta_Qs"]), 3)
        self.assertIn("-NQ_r*DeltaH/(X*H_r*H)", rows["Delta_Qs"])

    def test_both_pressure_terms_survive_Ns(self):
        rows = map_ledger()
        self.assertEqual(rows["Delta_Ns"][-2:], ("4*A*eta*DeltaCp", "-d*DeltaCp_eta"))

    def test_p2_denominator_difference_survives(self):
        rows = map_ledger()
        self.assertEqual(len(rows["Delta_p2"]), 2)
        self.assertIn("-X*Ns_r*DeltaE/(L*E_r*E)", rows["Delta_p2"])

    def test_R_power_table_covers_exact_moment_rows(self):
        rows = moment_ledger()
        envelope = envelope_ledger()
        self.assertLessEqual(max(row["R_power"] for row in rows), envelope["moments_over_t2_ea_s6"])
        self.assertEqual(envelope["Delta_p2_over_C_t_ea_s3"], 24)

    def test_worst_cone_factor_is_ts_times_p2(self):
        rows = cone_ledger()
        worst = max((row["C_power"], row["R_power"], row["name"])
                    for row in rows if row["C_power"] == 2)
        self.assertEqual(worst, (2, 27, "ts_Delta_p2"))
        self.assertTrue(all(row["C_power"] <= 2 and row["R_power"] <= 80 for row in rows))

    def test_selected_scale_inequalities_are_strict(self):
        data = scale_certificate()
        self.assertLess(data["log_R80_over_C_upper_per_B"], 0)
        self.assertLess(data["log_tR3_upper_per_B"], 0)
        self.assertLess(data["log_Kact_t_over_Cminus17_upper_per_B"], 0)
        self.assertEqual(data["cone_C_power_after_R160_below_C2"], -31)

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
