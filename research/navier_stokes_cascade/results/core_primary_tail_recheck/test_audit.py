#!/usr/bin/env python3
"""Exact coupled-gate, band-overlap, edge and evidence tests."""

from copy import deepcopy
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import (active_bands, band_checks, dyadic, edge_constant,
                          fixture, flat_constant, inherited, physical_losses,
                          review, select_tail)


class PrimaryTailRecheckTests(unittest.TestCase):
    def test_complete_actual_budget_keys_required(self):
        p = fixture()
        p.pop("length")
        with self.assertRaises(ValueError):
            select_tail(p)

    def test_floating_budgets_rejected(self):
        p = fixture()
        p["h"] = .005
        with self.assertRaises(TypeError):
            select_tail(p)

    def test_construction_ranges_checked(self):
        for key, value in (("h", F(0)), ("radius", F(1)), ("epsilon", F(2)), ("slow_cost", F(1, 2))):
            p = fixture()
            p[key] = value
            with self.assertRaises(ValueError):
                select_tail(p)

    def test_all_primary_gates_on_one_integer(self):
        selected = select_tail(fixture())
        self.assertTrue(all(selected["checks"].values()))
        self.assertEqual(selected["ell_tail"], max(selected["thresholds"].values()))

    def test_later_bands_remain_accepted(self):
        p = fixture()
        row = select_tail(p)
        for delta in (0, 1, 100, 1000000):
            self.assertTrue(all(band_checks(p, row["ell_tail"] + delta, row["gamma"], row["binary_exponent_required"], row["slow_fourth_power_required"]).values()))

    def test_mixed_gate_agrees_with_inherited_small_fixture(self):
        p = fixture()
        old = inherited.tail_certificate(*(p[k] for k in ("h", "cost", "rb", "omega", "lpsi", "cmesh", "cpert", "dbox", "dchart", "epsilon")))
        row = select_tail(p)
        self.assertEqual(row["gamma"], old["gamma"])
        self.assertEqual(row["pointwise_budget"], old["pointwise_budget"])
        self.assertGreaterEqual(row["ell_tail"], old["ell_tail"])

    def test_pulse_length_gate_is_independent(self):
        p = fixture()
        p["length"] = F(10**40)
        row = select_tail(p)
        bound = row["thresholds"]["pulse_length"]
        self.assertGreaterEqual(2 * p["radius"] * bound**2, p["length"])
        self.assertLess(2 * p["radius"] * (bound - 1)**2, p["length"])

    def test_padded_slow_gate_has_exact_fourth_root_threshold(self):
        row = select_tail(fixture())
        n = row["thresholds"]["padded_slow_coordinates"]
        self.assertGreaterEqual(F(n**4), row["slow_fourth_power_required"])
        self.assertLess(F((n - 1)**4), row["slow_fourth_power_required"])

    def test_wrong_preliminary_domain_Q_alone_is_not_enough(self):
        qstar = F(1, 10)
        Q = 3 * qstar / 4
        self.assertLess(Q, qstar)
        self.assertGreater(4 * Q, qstar)

    def test_domain_enlargement_stays_inside_background(self):
        for E in (16, 32, 64):
            qstar = F(1, 2**90)
            q = min(qstar / 64, dyadic(-E - 2))
            for ratio in (F(1, 4), F(1), F(4)):
                Q = q / ratio
                self.assertLessEqual(Q, dyadic(-E))
                self.assertLessEqual(4 * Q, qstar / 4)

    def test_enormous_scale_exponent_is_only_a_short_integer_certificate(self):
        row = select_tail(fixture(F(1, 2**200)))
        self.assertTrue(all(row["checks"].values()))
        self.assertGreater(row["binary_scale_exponent"], 10**20)
        self.assertLess(row["binary_scale_exponent"].bit_length(), 500)
        self.assertNotIn("scale_floor", row)

    def test_no_derivative_order_enters_the_selector(self):
        p = fixture()
        p["derivative_order"] = 3
        with self.assertRaises(ValueError):
            select_tail(p)

    def test_five_closed_enlarged_bands_at_dyadic_point(self):
        self.assertEqual(active_bands(dyadic(-20)), [18, 19, 20, 21, 22])

    def test_overlap_bound_away_from_dyadic_points(self):
        for q in (F(1, 3), F(17, 19), F(123, 10000), dyadic(-100) * F(3, 7)):
            self.assertLessEqual(len(active_bands(q)), 5)
            self.assertTrue(all(F(1, 4) <= q / dyadic(-ell) <= 4 for ell in active_bands(q)))

    def test_band_sum_has_no_scale_count_loss(self):
        for q in (dyadic(-20), F(3, 7) * dyadic(-20)):
            for K in (0, 1, 3, 11):
                value = sum((dyadic(-ell)**K for ell in active_bands(q)), F(0))
                self.assertLessEqual(value, 5 * 4**K * q**K)

    def test_envelope_without_label_support_separation_has_no_fixed_count_bound(self):
        # Six coincident unit-normalized errors at one band violate the
        # five-term K=0 counting bound despite identical individual bounds.
        self.assertGreater(sum([F(1)] * 6), 5)
        self.assertFalse(audit.FLAGS["all_generated_sources_satisfy_inverse_input_hypotheses"])

    def test_exact_square_completion_includes_peak(self):
        a = F(1, 101)
        row = flat_constant(a, 11, F(13, 7), 7)
        center = row["D"] / (2 * a)
        self.assertEqual(row["D"] * center - a * center**2, row["log_constant"])

    def test_negative_physical_loss_keeps_valid_nonnegative_D(self):
        self.assertEqual(flat_constant(F(1), 2, F(-100), 7)["D"], 2)

    def test_factorial_majorant_keeps_each_edge_exponential(self):
        for j in range(7):
            for x in (F(1, 100), F(1), F(17)):
                exp_lower = sum((x**i / factorial(i) for i in range(j + 1)), F(0))
                self.assertLessEqual(x**j, factorial(j) * exp_lower)

    def test_original_narrow_edge_constant_cannot_be_replaced_by_one(self):
        self.assertGreater(edge_constant(8, F(1, 20000), F(2)), edge_constant(8, F(1), F(2)))
        self.assertEqual(edge_constant(8, F(1, 20000), F(2)), 3 + factorial(4) * (40000**4 + 1))

    def test_physical_phase_partitions_and_coefficient_costs(self):
        h, kappa = F(1, 200), F(1, 100000)
        for a in range(8):
            for b in range(8 - a):
                row = physical_losses(h, a, b)
                for j in range(a + b + 1):
                    phase = F(a, 2) + b * (1 + h) + j * h / 2
                    self.assertLessEqual(phase, a * row["sx"] + b * row["st"])
                self.assertLessEqual(F(1, 2) + h * kappa, row["sx"])
                self.assertLessEqual(1 + h, row["st"])

    def test_every_review_and_inherited_source_gate(self):
        self.assertTrue(all(review()["checks"].values()))
        record = audit.generate()
        self.assertTrue(all(record["gates"].values()))
        self.assertGreater(len(record["inherited_source_checks"]), 80)

    def test_false_scope_and_frontier_changes_rejected(self):
        record = audit.generate()
        for key, value in audit.FLAGS.items():
            if not value:
                changed = deepcopy(record)
                changed["flags"][key] = True
                with self.assertRaises(ValueError):
                    audit.validate(changed)
        changed = deepcopy(record)
        changed["remaining_obligations"].pop()
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_frozen_record_and_fresh_replay(self):
        audit.validate(json.loads((audit.HERE / "evidence.json").read_text()))
        with tempfile.TemporaryDirectory(dir=audit.PROJECT) as directory:
            path = Path(directory) / "fresh.json"
            audit.main(["--verify-record", str(audit.HERE / "evidence.json"), "--output", str(path)])
            audit.validate(json.loads(path.read_text()))
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])


if __name__ == "__main__":
    unittest.main()
