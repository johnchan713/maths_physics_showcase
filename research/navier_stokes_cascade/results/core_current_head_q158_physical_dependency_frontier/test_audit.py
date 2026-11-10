import json
import tempfile
import unittest
from pathlib import Path
import audit
import q158


class Q158Tests(unittest.TestCase):
    def test_route_is_current_and_ordered(self):
        route = q158.physical_route()
        self.assertEqual(len(route), 10)
        self.assertIn("Q147", route[0])
        self.assertIn("Q157", route[-1])

    def test_complete_primitive_frontier(self):
        rows = q158.primitive_frontier()
        self.assertEqual([row["id"] for row in rows], [f"P{i}" for i in range(1, 8)])
        self.assertEqual(rows[0]["status"], "discharged-by-Q158-provider-binding")
        self.assertTrue(all(row["status"].startswith("conditional") for row in rows[1:]))

    def test_first_binding(self):
        row = q158.first_binding()
        self.assertTrue(row["same_profile_required"])
        self.assertTrue(row["discharged"])
        self.assertFalse(row["numerical_materialization_required"])
        self.assertIn("V0<=N", row["conclusion"])

    def test_quantifier_order(self):
        rows = q158.quantifier_order()
        self.assertEqual(len(rows), 8)
        self.assertIn("then choose that Borel scale", rows[1])
        self.assertIn("T<1", rows[-1])
        self.assertIn("R to infinity", rows[-1])

    def test_provider_contracts(self):
        row = q158.provider_contracts()
        self.assertTrue(all(row.values()))

    def test_record_scope(self):
        row = q158.record()
        self.assertTrue(row["first_loose_provider_binding_discharged"])
        self.assertFalse(row["unconditional_admissible_force_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_identity_guard(self):
        for fn in (q158.physical_route, q158.primitive_frontier,
                   q158.first_binding, q158.quantifier_order,
                   q158.provider_contracts):
            with self.assertRaises(ValueError):
                fn("wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 395)
        self.assertGreater(record["closure"]["edges"], 1270)
        self.assertGreater(record["closure"]["source_bindings"], 5365)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory() as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
