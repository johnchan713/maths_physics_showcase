import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q155


class Q155Tests(unittest.TestCase):
    def test_exact_twelve_term_ledger(self):
        ledger = q155.localization_ledger()
        self.assertEqual(ledger["count"], 12)
        self.assertEqual(len(set(ledger["terms"])), 12)
        self.assertEqual(ledger["viscosity"], 1)
        self.assertEqual(ledger["complete_pressure"], "p0 grad c")
        self.assertTrue(ledger["force_divergence_retained"])

    def test_thirteenth_term_vanishes(self):
        self.assertIn("grad(c) cross A", q155.localization_ledger()["vanishing_thirteenth"])

    def test_derivative_ledger(self):
        r = q155.derivative_requests(7)
        self.assertEqual(r["force"], {"potential": 9, "direct_swirl": 8, "pressure": 7})
        self.assertEqual(r["cauchy"], {"potential": 10, "direct_swirl": 9,
                                        "pressure": 8, "residual": 8})
        self.assertEqual(r["Q154_spare"], {"potential": 1, "direct_swirl": 1,
                                            "pressure": 1})

    def test_corner_transfer(self):
        r = q155.corner_transfer(5)
        self.assertEqual(r["rho_power"], 6)
        self.assertEqual(r["constant_factor"], 4 ** 6)
        self.assertTrue(r["cartesian_flat_for_arbitrary_N"])

    def test_borel_power_decays_on_admitted_diagonal(self):
        for j in range(2, 30):
            for a in range(j // 2 + 1):
                for m in range(j // 2 - a + 1):
                    if q155.borel_admitted(j, a, m):
                        self.assertLess(q155.borel_term_power(j, m), 0)

    def test_energy_budget(self):
        b = q155.energy_budget(Q(7, 3))
        self.assertEqual(b["sup_energy"], Q(49, 9))
        self.assertEqual(b["integrated_dissipation"], Q(49, 18))
        self.assertFalse(b["sum_of_separate_extrema_budget_claimed"])

    def test_provider_contracts(self):
        c = q155.provider_contracts()
        self.assertTrue(c["Q154_supplies_Q8_punctured_endpoint_premise"])
        self.assertTrue(c["Q153_supplies_corner_complete_residual_flatness"])
        self.assertFalse(c["future_velocity_or_pressure"])

    def test_record_scope(self):
        r = q155.record()
        self.assertTrue(r["compact_force_and_energy_certified_conditionally"])
        self.assertFalse(r["future_solution_pressure_comparison_or_blowup_certified"])

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            q155.localization_ledger("wrong")
        with self.assertRaises(ValueError):
            q155.derivative_requests(-1)
        with self.assertRaises(ValueError):
            q155.corner_transfer(-1)
        with self.assertRaises(ValueError):
            q155.energy_budget(-1)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 390)
        self.assertGreater(record["closure"]["source_bindings"], 5300)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory() as td:
            target = Path(td) / "fresh.json"
            data = audit.generate()
            target.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), data)


if __name__ == "__main__":
    unittest.main()
