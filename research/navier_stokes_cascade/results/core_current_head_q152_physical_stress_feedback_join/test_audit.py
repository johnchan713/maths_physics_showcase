import json
import tempfile
import unittest
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import audit
import q152 as q


class Tests(unittest.TestCase):
    def test_physical_scale_identity(self):
        out = q.scale_identity()
        self.assertEqual(out["A"], Q(101, 200))
        self.assertEqual(out["band_power"], 0)
        self.assertTrue(out["exact"])

    def test_same_family_root_and_squared_partition(self):
        # Each label has (eta^2,a_scale,H,y,target).  Both root equations
        # and the physical normalization are exact in this rational fixture.
        labels = [
            (Q(1, 4), 2, 3, 5, 15),
            (Q(3, 4), 2, 3, 5, 15),
        ]
        out = q.principal_assembly(labels, 30)
        self.assertTrue(out["exact"])
        self.assertEqual(out["cover"], 1)
        self.assertEqual(out["defect"], 0)

    def test_root_scale_and_cover_errors_are_detected(self):
        good = [(Q(1, 4), 2, 3, 5, 15), (Q(3, 4), 2, 3, 5, 15)]
        bad_root = deepcopy(good); bad_root[0] = (Q(1, 4), 2, 4, 5, 15)
        bad_cover = good[:1]
        self.assertFalse(q.principal_assembly(bad_root, 30)["exact"])
        self.assertFalse(q.principal_assembly(bad_cover, 30)["exact"])

    def test_complete_curl_bound_has_half_power_margin(self):
        out = q.rates()
        self.assertEqual(out["required_stress_after_epsilon_inverse"], Q(24999, 50000))
        self.assertEqual(out["available_complete_curl_after_epsilon_inverse"], Q(49999, 50000))
        self.assertEqual(out["complete_curl_margin"], Q(1, 2))
        self.assertEqual(out["feedback"], out["D"])

    def test_q151_gate_is_closed_without_family_identification(self):
        out = q.stress_gate()
        self.assertTrue(out["principal_gate_closed_exactly"])
        self.assertFalse(out["polynomial_family_closeness_used"])
        self.assertEqual(out["complete_curl_margin"], Q(1, 2))

    def test_physical_provider_join(self):
        out = q.provider_contracts()
        self.assertTrue(out["same_physical_family"])
        self.assertFalse(out["canonical_and_physical_entrances_identified"])
        self.assertEqual(out["mean_feedback"], "D")
        self.assertIn("pole 9p", out["next_regular_source"])
        self.assertIn("A1", out["fresh_fast"])

    def test_missing_physical_contract_rejected(self):
        original = q.evidence
        try:
            records = {name: original(name) for name in q.PROVIDERS}
            records["actual_band_pulse_bounds"] = deepcopy(records["actual_band_pulse_bounds"])
            records["actual_band_pulse_bounds"]["flags"]["conditional_real_primary_curl_covariance_remainder_gain_supplied"] = False
            q.evidence = lambda name: records[name]
            with self.assertRaises(ValueError):
                q.provider_contracts()
        finally:
            q.evidence = original

    def test_invalid_inputs_rejected(self):
        for k in (0, Q(1, 4), -1):
            with self.assertRaises(ValueError):
                q.rates(k)
        for h in (0, Q(1, 100), -1):
            with self.assertRaises(ValueError):
                q.scale_identity(h)
        with self.assertRaises(ValueError):
            q.principal_assembly([], 0)
        with self.assertRaises(ValueError):
            q.provider_contracts("other")

    def test_scope_record(self):
        record = q.record()
        self.assertFalse(record["canonical_and_physical_primary_families_identified"])
        self.assertFalse(record["generated_source_endpoint_and_complete_tail_join_certified"])
        self.assertFalse(record["complete_residual_or_blowup_certified"])
        self.assertIn("endpoint extension", record["next"])

    def test_frozen_record(self):
        with tempfile.TemporaryDirectory() as directory:
            fresh = Path(directory) / "fresh.json"
            fresh.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(fresh.read_bytes(), (audit.HERE / "evidence.json").read_bytes())


if __name__ == "__main__":
    unittest.main()
