from copy import deepcopy
from fractions import Fraction as Q
import json, unittest
from pathlib import Path
import audit
from closure import (closure_statement, fixed_order_induction,
                     normalized_request, physical_interface,
                     physical_request, radial_rate_closure, topological_order)


class InitialA0PhysicalTailClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def control(self, key):
        found = [value for name, value in self.record["review"]["checks"].items() if key in name]
        self.assertTrue(found, key)
        self.assertTrue(all(found), key)

    def test_acyclic_and_noncircular(self):
        self.control("acyclic")
        self.control("not_a_gr_input")
        self.assertEqual(topological_order()[-1], "A0_Qcirc_gr_final")

    def test_normalized_requests(self):
        self.control("source_and_graph_indices")
        self.control("latest_A1_route")
        for n in (0, 1, 9):
            r = normalized_request(n)
            self.assertEqual(r["A0_source_order"], 399995+n)
            self.assertEqual(r["largest_same_prefix_order"], 400009+n)

    def test_normalized_rates(self):
        self.control("strict_kappa_margin")
        self.assertEqual(radial_rate_closure()["strict_margin"], Q(1,100000))

    def test_physical_request_and_interface(self):
        self.control("physical_zero_jet")
        self.control("positive_margin")
        self.assertGreater(physical_request()["A0_strict_margin"], 0)
        self.assertFalse(physical_interface()["complete_regular_residual_improvement_proved"])

    def test_finite_order_induction(self):
        self.control("covers_normalized_and_physical")
        self.assertEqual(fixed_order_induction(physical=True)["request"]["largest_same_prefix_order"], 60200021)

    def test_scope(self):
        self.control("conditionally_not_unconditionally")
        self.control("scope_is_honest")
        self.assertFalse(closure_statement()["complete_residual_or_blowup_closed"])

    def test_guards(self):
        for bad in (0, -Q(1,1000), Q(1,50)):
            with self.assertRaises(ValueError):
                normalized_request(kappa=bad)
        with self.assertRaises(ValueError):
            physical_request(h=Q(1,50))

    def test_frozen_record_and_tampering(self):
        audit.validate(json.loads((Path(__file__).parent/"evidence.json").read_text()))
        for key in ("flags", "review", "source_hashes"):
            record = deepcopy(self.record)
            record[key]["tamper"] = True
            with self.assertRaises(ValueError):
                audit.validate(record)


if __name__ == "__main__":
    unittest.main()
