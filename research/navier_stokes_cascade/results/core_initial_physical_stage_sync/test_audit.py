from copy import deepcopy
from fractions import Fraction as Q
import json, unittest
from pathlib import Path
import audit
from stage_sync import compact_difference_identity, manufactured_change, scope, stages


class InitialPhysicalStageSyncTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_exact_stage_orders_and_exponents(self):
        s = stages()
        self.assertEqual(s["first_real"]["regular_order"], Q(99999, 100000))
        self.assertEqual(s["first_sourced"]["regular_order"], Q(74999, 50000))
        self.assertEqual(s["first_real"]["physical_exponent"], Q(-30100001, 20000000))
        self.assertEqual(s["first_sourced"]["physical_exponent"], Q(-15025001, 10000000))
        self.assertGreater(s["physical_gain"], 0)

    def test_neither_synchronized_row_has_absolute_decay(self):
        for m in (0, 1, 4):
            s = stages(m=m, N=0)
            self.assertLess(s["first_real"]["physical_exponent"], 0)
            self.assertLess(s["first_sourced"]["physical_exponent"], 0)
            self.assertFalse(s["first_sourced_meets_requested_bound"])

    def test_compact_state_changes_exactly(self):
        identity = compact_difference_identity()
        self.assertIn("Delta_gr", identity["difference"])
        self.assertIn("Delta_Phi", identity["difference"])
        self.assertFalse(identity["initial_compact_bound_transfers_to_sourced_state"])

    def test_manufactured_noncancellation(self):
        f = manufactured_change()
        self.assertEqual(f["Qcirc_gr"], (Q(3, 2), Q(-3, 2)))
        self.assertEqual(f["Qcirc_phi"], (Q(7), Q(-7)))
        self.assertTrue(f["radial_changed"] and f["axial_changed"])

    def test_scope(self):
        flags = scope()
        self.assertFalse(flags["synchronized_post_sourced_compact_DAG_closed_here"])
        self.assertFalse(flags["complete_physical_residual_decay_proved"])
        self.assertFalse(flags["blowup_proved"])

    def test_guards(self):
        for bad_h in (0, Q(1, 100), -Q(1, 10)):
            with self.assertRaises(ValueError):
                stages(h=bad_h)
        with self.assertRaises(ValueError):
            stages(m=-1)
        with self.assertRaises(ValueError):
            stages(N=-1)

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
