import json
import tempfile
import unittest
from pathlib import Path
import audit
import q159


class Q159Tests(unittest.TestCase):
    def test_P2_is_same_object_and_closed(self):
        row = q159.p2_contract()
        self.assertTrue(row["analytically_discharged"])
        self.assertTrue(row["complete_pressure_and_cutoff_curl_retained"])
        self.assertFalse(row["new_profile_or_scale_sequence_selected"])
        self.assertIn("before", row["coefficient_quantifier"])
        self.assertIn("derivative-independent", row["pulse_quantifier"])

    def test_P3_gate_is_universal_and_open(self):
        row = q159.p3_gate()
        self.assertTrue(row["obstruction_open"])
        self.assertTrue(row["finite_stage_result_is_conditional_on_this_request"])
        self.assertFalse(row["source_envelope_provider_marks_every_later_target_verified"])
        self.assertIn("every finite cycle", row["universal_request"])
        self.assertIn("zeta", row["weighted_bound"])

    def test_pressure_and_moments_are_complete(self):
        row = q159.p3_gate()
        self.assertIn("J_theta", row["angular_moment"])
        self.assertIn("c_rho P_mean", row["axial_moment"])
        for required in ("axial velocity", "six covariance", "Pi1/Pi2", "zero-Haar"):
            self.assertIn(required, row["memory"])

    def test_corrected_frontier(self):
        rows = q159.corrected_frontier()
        self.assertEqual([row["id"] for row in rows], [f"P{i}" for i in range(1, 8)])
        self.assertEqual(rows[1]["status"], "closed-by-Q149-and-Q159-audit")
        self.assertEqual(rows[2]["status"], "open-at-all-finite-cycle-target-quantifier")

    def test_provider_contracts(self):
        self.assertTrue(all(q159.provider_contracts().values()))

    def test_record_scope(self):
        row = q159.record()
        self.assertTrue(row["P2_analytically_discharged"])
        self.assertFalse(row["P3_universal_finite_cycle_gate_discharged"])
        self.assertFalse(row["unconditional_admissible_force_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_identity_guard(self):
        for fn in (q159.p2_contract, q159.p3_gate,
                   q159.corrected_frontier, q159.provider_contracts):
            with self.assertRaises(ValueError):
                fn("wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 398)
        self.assertGreater(record["closure"]["edges"], 1280)
        self.assertGreater(record["closure"]["source_bindings"], 5410)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory() as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
