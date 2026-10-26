from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import unittest
import audit
from recurrence import (closure_statement, first_stage_for_power, layer,
                        physical_request, pole_support_ledger, rates,
                        state_ceiling)


class FiniteStageCompactRecurrenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_rates(self):
        self.assertEqual(rates(0)["regular"], Q(74999, 50000))
        self.assertEqual(rates(2)["regular"]-rates(1)["regular"], Q(49997, 100000))

    def test_actual_layer_offsets(self):
        r = layer(0)
        self.assertEqual(r["sourced_backward_offset"], 200018)
        self.assertEqual(r["generic_backward_offset"], 300011)
        self.assertEqual(r["full_cycle_backward_offset"], 500029)

    def test_state_recurrence(self):
        self.assertEqual(state_ceiling(0)["first_real_ceiling"], 500038)
        self.assertEqual(state_ceiling(1)["first_real_ceiling"], 1000067)
        self.assertGreater(state_ceiling(3)["first_real_ceiling"], state_ceiling(2)["first_real_ceiling"])

    def test_physical_seed_and_next_stage(self):
        self.assertEqual(physical_request(0)["first_real_ceiling"], 60500042)
        r = physical_request(1)
        self.assertEqual(r["first_real_ceiling"], 60900077)
        self.assertEqual(min(r["strict_margins"].values()), Q(1, 20000000))

    def test_fixed_stage_poles_and_supports(self):
        r = pole_support_ledger(3)
        self.assertTrue(r["finite_for_fixed_stage"])
        self.assertEqual(r["harmonics_at_most"], 32)
        self.assertEqual(len(r["distinct_common_poles"]), 5)
        self.assertFalse(r["actual_constants_or_high_jets_materialized"])

    def test_physical_threshold(self):
        r = first_stage_for_power()
        self.assertEqual(r["stage"], 1002)
        self.assertGreaterEqual(r["achieved_power"], 0)
        self.assertLess(r["preceding_power"], 0)

    def test_guards(self):
        with self.assertRaises(ValueError):
            rates(-1)
        with self.assertRaises(ValueError):
            physical_request(0, N=-1)
        with self.assertRaises(ValueError):
            pole_support_ledger(0, original_pole=-1)
        with self.assertRaises(ValueError):
            first_stage_for_power(h=Q(1, 100))

    def test_scope(self):
        r = closure_statement()
        self.assertTrue(r["finite_stage_compact_graph_synchronized_conditionally"])
        self.assertFalse(r["one_field_all_targets_diagonal_proved"])
        self.assertFalse(r["localization_forcing_energy_or_blowup_proved"])

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
