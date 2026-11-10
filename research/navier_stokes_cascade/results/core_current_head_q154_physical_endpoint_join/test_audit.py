import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q154


class Q154Tests(unittest.TestCase):
    def test_selected_parameters(self):
        p = q154.selected_parameters()
        self.assertEqual(p["D"], Q(99, 200))
        self.assertEqual(p["implicit_floor"], Q(99, 100))
        self.assertGreater(p["stage_gain"], 0)

    def test_full_enlarged_band_window(self):
        w = q154.band_window(24, 20, 1)
        self.assertEqual((w["lower"], w["upper"]), (18, 26))
        self.assertEqual(len(w["bands"]), 9)

    def test_narrow_window_would_omit_edges(self):
        w = q154.band_window(24, 20, 1)
        self.assertIn(18, w["bands"])
        self.assertIn(26, w["bands"])

    def test_observation_margin(self):
        m = q154.observation_margin(Q(1, 64), 10)
        self.assertTrue(m["inside_qstar"])
        self.assertTrue(m["first_band_margin"])
        self.assertEqual(m["sixteen_qmax"], 6 * m["qe"])

    def test_endpoint_requests(self):
        r = q154.endpoint_requests(3, 5)
        self.assertEqual(r, {"m": 8, "potential": 12, "direct_swirl": 11,
                             "pressure": 10, "corner_residual": 9})

    def test_complete_pressure_convolution(self):
        h = Q(1, 200)
        e = q154.exterior_coefficients(2, h)
        A1 = q154.exterior_coefficients(1, h)["A_j"]
        A2 = e["A_j"]
        self.assertEqual(e["S_j"] - 2 * A2, 2 * A1 * A1)

    def test_taylor_ratio_is_unbounded(self):
        vals = [q154.taylor_ratio(j) for j in (10, 100, 1000)]
        self.assertLess(vals[0], vals[1])
        self.assertLess(vals[1], vals[2])
        self.assertGreater(vals[2], 4000)

    def test_provider_contracts(self):
        c = q154.provider_contracts()
        self.assertTrue(c["q7_blocker_supplied_conditionally_by_later_physical_chain"])
        self.assertTrue(c["same_Q153_field_cutoffs_supports_and_pressure_gauge"])
        self.assertFalse(c["singular_origin_included"])

    def test_record_scope(self):
        r = q154.record()
        self.assertTrue(r["punctured_closed_time_endpoint_certified_conditionally"])
        self.assertFalse(r["singular_origin_regularized"])
        self.assertFalse(r["compact_force_energy_comparison_or_blowup_certified"])

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            q154.selected_parameters(h=0)
        with self.assertRaises(ValueError):
            q154.band_window(10, 20, 1)
        with self.assertRaises(ValueError):
            q154.endpoint_requests(-1, 0)
        with self.assertRaises(ValueError):
            q154.taylor_ratio(0, r=0)

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
