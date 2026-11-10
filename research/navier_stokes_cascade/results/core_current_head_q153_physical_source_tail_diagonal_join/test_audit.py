import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q153


class Q153Tests(unittest.TestCase):
    def test_rates(self):
        r = q153.rates()
        self.assertEqual(r["D"], Q(3, 2) - Q(2, 100000))
        self.assertEqual(r["delta"], Q(1, 2) - Q(3, 100000))
        self.assertGreater(r["delta"], 0)

    def test_bell_quotient(self):
        out = q153.bell_quotient([2, 3, 5, 7], [0, 11, 13, 17])
        self.assertEqual(out["inverse_gaussian"][:3], [Q(1), Q(11), Q(134)])
        self.assertEqual(out["quotient"][0], 2)
        self.assertEqual(out["quotient"][1], 25)
        self.assertEqual(out["quotient"][2], 339)

    def test_quotient_rejects_negative(self):
        with self.assertRaises(ValueError):
            q153.bell_quotient([1, -1], [0, 1])

    def test_physical_request_is_strict_and_finite(self):
        req = q153.compact_request(0, 0)
        self.assertEqual(req["p0"], 30100002)
        self.assertEqual(req["p1"], 30050003)
        self.assertEqual(req["J0"], 60200007)
        self.assertEqual(req["J1"], 60100009)
        self.assertEqual(req["first_real_ceiling"], 60400023)
        self.assertGreater(req["radial_margin"], 0)
        self.assertGreater(req["axial_margin"], 0)

    def test_requests_increase(self):
        a = q153.compact_request(0, 0)
        b = q153.compact_request(2, 3)
        self.assertGreater(b["J0"], a["J0"])
        self.assertGreater(b["first_real_ceiling"], a["first_real_ceiling"])

    def test_cutoff_allocation(self):
        c = q153.cutoff_allocation()
        self.assertTrue(c["double_count_free"])
        self.assertIn("pressure gradients", c["per_stage_regular"])
        self.assertIn("A0/A1", c["per_stage_additive"])
        self.assertIn("scale cutoffs", c["diagonal"])

    def test_provider_contracts(self):
        c = q153.provider_contracts()
        self.assertIn("Q152", c["physical_seed"])
        self.assertFalse(c["independent_continuum_validation"])

    def test_record_scope(self):
        r = q153.record()
        self.assertTrue(r["physical_generated_source_endpoint_join_certified_conditionally"])
        self.assertTrue(r["all_fixed_stage_additive_physical_tails_certified_conditionally"])
        self.assertFalse(r["closed_time_forcing_energy_or_blowup_certified"])

    def test_wrong_identity_rejected(self):
        with self.assertRaises(ValueError):
            q153.rates(identity="wrong")
        with self.assertRaises(ValueError):
            q153.cutoff_allocation(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 300)
        self.assertGreater(record["closure"]["source_bindings"], 4000)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory() as td:
            target = Path(td) / "fresh.json"
            data = audit.generate()
            target.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), data)


if __name__ == "__main__":
    unittest.main()
