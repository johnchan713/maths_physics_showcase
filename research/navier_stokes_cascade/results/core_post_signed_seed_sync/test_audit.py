from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import unittest
import audit
from ledger import (corrected_seed, interface, normalized_seed_request,
                    physical_seed_request, stage_obligation)


class PostSignedSeedSyncTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_normalized_expansion(self):
        r = normalized_seed_request()
        self.assertEqual(r["first_sourced_target_order"], 200017)
        self.assertEqual(r["first_real_ceiling"], 500038)

    def test_physical_expansion(self):
        r = physical_seed_request()
        self.assertEqual(r["first_sourced_target_order"], 60200021)
        self.assertEqual(r["first_real_ceiling"], 60500042)
        self.assertEqual(Q(r["illustrative_domain"]["q"]), 361200138)

    def test_physical_margins(self):
        r = physical_seed_request()
        self.assertEqual(r["current_A0"]["strict_margin"], Q(1, 20000000))
        self.assertEqual(r["current_A1"]["strict_margin"], Q(1, 20000000))

    def test_corrected_seed(self):
        r = corrected_seed()
        self.assertEqual(r["seed_margin_over_old_minimum"], Q(1, 25000))
        self.assertEqual(r["mean_relation"], Q(1, 2))
        self.assertTrue(r["half_line_rate_recurrence_applicable"])

    def test_later_stage_obligation_is_explicit(self):
        r = stage_obligation()
        self.assertGreater(r["strict_margins"]["A0"], 0)
        self.assertGreater(r["strict_margins"]["A1"], 0)
        self.assertEqual(len(r["required_inequalities"]), 2)
        self.assertFalse(r["stage_indexed_backward_request_recurrence_supplied"])

    def test_scope(self):
        r = interface()
        self.assertFalse(r["old_half_line_rate_recurrence_alone_proves_later_compact_inputs"])
        self.assertFalse(r["complete_finite_stage_residual_decay_proved"])
        self.assertFalse(r["infinite_realization_or_blowup_proved"])

    def test_guards(self):
        with self.assertRaises(ValueError):
            normalized_seed_request(-1)
        with self.assertRaises(ValueError):
            physical_seed_request(N=-1)
        with self.assertRaises(ValueError):
            stage_obligation(0)
        with self.assertRaises(ValueError):
            stage_obligation(h=Q(1, 100))

    def test_frozen_record_and_tampering(self):
        path = Path(__file__).parent/"evidence.json"
        if path.exists():
            audit.validate(json.loads(path.read_text()))
        for key in ("flags", "review", "source_hashes"):
            record = deepcopy(self.record)
            record[key]["tamper"] = True
            with self.assertRaises(ValueError):
                audit.validate(record)


if __name__ == "__main__":
    unittest.main()
