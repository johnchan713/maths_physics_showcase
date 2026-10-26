from copy import deepcopy
from fractions import Fraction as Q
import json, unittest
from pathlib import Path
import audit
from closure import (fixed_order_induction, normalized_request,
                     physical_request, residual_interface,
                     topological_order)


class FirstSourcedPhysicalTailClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_graph_order(self):
        order = topological_order()
        self.assertLess(order.index("postwave_A1_response"), order.index("fresh_fast_means"))
        self.assertEqual(order[-1], "A0_Qcirc_gr_S")

    def test_normalized_safe_request(self):
        r = normalized_request()
        self.assertEqual(r["raw_sources"], {"A0": 100005, "A1": 7})
        self.assertEqual(r["largest_first_real_ceiling"], 300021)

    def test_physical_safe_request(self):
        r = physical_request()
        self.assertEqual(r["raw_sources"], {"A0": 60200007, "A1": 60100009})
        self.assertEqual(r["direct_A1_first_real_ceiling"], 60300013)
        self.assertEqual(r["radial_nested_A1_first_real_ceiling"], 60400023)
        self.assertTrue(r["radial_branch_dominates"])

    def test_compact_strict_margins(self):
        r = physical_request()
        self.assertEqual(r["A0_physical_exponent"], Q(1, 20000000))
        self.assertEqual(r["A1_physical_exponent"], Q(1, 20000000))
        self.assertTrue(residual_interface()["compact_flat_for_each_fixed_m_N"])

    def test_regular_obstruction_remains(self):
        r = residual_interface()
        self.assertEqual(r["regular_physical_exponent"], Q(-15025001, 10000000))
        self.assertFalse(r["current_regular_decay_certified"])
        self.assertFalse(r["complete_residual_decay_proved"])

    def test_fixed_order_induction_scope(self):
        r = fixed_order_induction(physical=True)
        self.assertTrue(r["requested_ceiling_is_finite"])
        self.assertFalse(r["desired_A0_output_used_earlier"])
        self.assertFalse(r["primitive_contracts_independently_certified"])

    def test_guards(self):
        with self.assertRaises(ValueError):
            physical_request(m=-1)
        with self.assertRaises(ValueError):
            physical_request(N=-1)
        with self.assertRaises(ValueError):
            physical_request(h=Q(1, 100))
        with self.assertRaises(ValueError):
            residual_interface(eta=0)

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
