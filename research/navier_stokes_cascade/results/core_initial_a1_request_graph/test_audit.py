from copy import deepcopy
from fractions import Fraction as Q
import json, unittest
from pathlib import Path
import audit
from requests import pole_domains, qualitative_closure, request_graph


class InitialA1RequestGraphTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def control(self, key):
        found = [value for name, value in self.record["review"]["checks"].items() if key in name]
        self.assertTrue(found, key)
        self.assertTrue(all(found), key)

    def test_strict_count_and_raw_orders(self):
        self.control("least_strict")
        self.control("DzPhi_and_Phi")
        self.control("strict_kappa")

    def test_complete_finite_padding(self):
        self.control("same_state_pressure")
        self.control("target_profile_and_Haar")
        self.control("extra_sharp_initial")
        self.control("canonical_covariance")

    def test_scope_and_blocker(self):
        self.control("does_not_import")
        self.control("profile_and_support")
        self.control("does_not_promote")
        self.control("precise_Phi")

    def test_monotone_output_shift(self):
        a, b = request_graph(0), request_graph(7)
        for key in ("raw_Dz_Phi_order", "raw_Phi_order", "conservative_same_state_order",
                    "signed_target_order", "fixed_profile_order", "largest_immediate_order"):
            self.assertEqual(b[key]-a[key], 7)

    def test_guards(self):
        for bad in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                request_graph(bad)
        for bad in (0, -Q(1, 1000), Q(1, 50)):
            with self.assertRaises(ValueError):
                request_graph(0, bad)
        with self.assertRaises(ValueError):
            pole_domains(0, old_pole=-1)

    def test_fixed_values(self):
        graph = request_graph()
        self.assertEqual(graph["largest_immediate_order"], 200018)
        self.assertEqual(pole_domains()["q"], 1200030)
        self.assertFalse(qualitative_closure()["qualitative_smoothness_implies_weighted_field_rate"])

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
