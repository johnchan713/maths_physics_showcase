from copy import deepcopy
from fractions import Fraction as Q
import json, unittest
from pathlib import Path
import audit
from closure import closure_statement, fixed_order_induction, rate_closure, topological_order


class InitialPhiHighJetClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def control(self, key):
        found = [value for name, value in self.record["review"]["checks"].items() if key in name]
        self.assertTrue(found, key)
        self.assertTrue(all(found), key)

    def test_acyclic_and_noncircular(self):
        self.control("acyclic")
        self.control("not_a_Phi_input")
        self.assertEqual(topological_order()[-1], "Qcirc_Phi")

    def test_finite_high_order_induction(self):
        self.control("200018")
        self.control("numerically_materialized")
        for n in (0, 1, 9):
            self.assertEqual(fixed_order_induction(n)["request"]["largest_immediate_order"], 200018+n)

    def test_rates(self):
        self.control("order_one")
        self.control("three_minus_five")
        self.control("strict_kappa")
        self.assertEqual(rate_closure()["strict_margin"], Q(1,100000))

    def test_scope_and_conditional_inequality(self):
        self.control("follow_conditionally")
        self.control("not_silently_promoted")
        self.control("scope_is_honest")
        self.assertFalse(closure_statement()["actual_primitive_contracts_certified_here"])

    def test_guards(self):
        for bad in (0, -Q(1,1000), Q(1,50)):
            with self.assertRaises(ValueError):
                rate_closure(bad)

    def test_frozen_record(self):
        audit.validate(json.loads((Path(__file__).parent/"evidence.json").read_text()))
        self.assertTrue(all(self.record["gates"].values()))

    def test_tampering(self):
        for key in ("flags", "review", "source_hashes"):
            record = deepcopy(self.record)
            record[key]["tamper"] = True
            with self.assertRaises(ValueError):
                audit.validate(record)


if __name__ == "__main__":
    unittest.main()
