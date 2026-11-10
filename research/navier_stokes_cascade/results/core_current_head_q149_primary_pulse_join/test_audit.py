import json
import tempfile
import unittest
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import q149 as q
import audit


class Tests(unittest.TestCase):
    def q148(self):
        return json.loads((q.HERE.parent / "core_current_head_q148_coefficient_support_pressure/evidence.json").read_text())["symbolic"]

    def test_rule_schedule_disposes_of_external_numeric_scale_blocker(self):
        row = q.schedule_row(7, 4, Q(7, 200), 123, 4000, q.IDENTITY)
        self.assertTrue(row["strict_padded_rule_holds"])
        self.assertFalse(row["preexisting_numeric_scale_certificate_required"])
        bad = q.schedule_row(7, 4, Q(7, 200), 123, 100, q.IDENTITY)
        self.assertFalse(bad["strict_padded_rule_holds"])

    def test_same_profile_provider_contracts_close(self):
        out = q.close(self.q148(), q.IDENTITY)
        self.assertTrue(out["uniform_primary_pulse_family_certified"])
        self.assertFalse(out["canonical_and_physical_primary_families_identified"])
        self.assertFalse(out["complete_nonlinear_generated_pulse_hierarchy_certified"])
        self.assertIn("Pi1", out["full_pressure_memory"])
        self.assertEqual(set(out["supports"]), {"transverse", "pulse", "radial", "slow"})

    def test_missing_same_profile_or_provider_contract_rejected(self):
        bad = self.q148()
        bad["status"] = "stage-only"
        with self.assertRaises(ValueError):
            q.close(bad, q.IDENTITY)
        original = q.provider_flags
        try:
            flags = original()
            flags = deepcopy(flags)
            flags["physical_pulse_transfer"]["actual_background_inputs_to_primary_pulse_theorem_discharged"] = False
            q.provider_flags = lambda: flags
            with self.assertRaises(ValueError):
                q.close(self.q148(), q.IDENTITY)
        finally:
            q.provider_flags = original

    def test_exact_physical_derivative_powers(self):
        self.assertEqual(q.primary_powers(Q(1, 200), 0, 0, q.IDENTITY), {
            "velocity": Q(-201, 400), "pressure": Q(-201, 200), "potential": Q(0)
        })
        self.assertEqual(q.primary_powers(Q(1, 200), 2, 1, q.IDENTITY), {
            "velocity": Q(-503, 200), "pressure": Q(-1207, 400), "potential": Q(-161, 80)
        })

    def test_covariance_obstruction_is_not_absolute_decay(self):
        out = q.covariance_obstruction(Q(1, 200), q.IDENTITY)
        self.assertEqual(out["Q_exponent"], Q(-10000001, 10000000))
        self.assertFalse(out["absolutely_decaying"])
        self.assertIn("axial velocity", out["necessary_next_step"])

    def test_scope_record(self):
        record = q.record()
        self.assertFalse(record["complete_covariance_residual_or_blowup_certified"])
        self.assertIn("primary family only", record["scope"])

    def test_frozen_record(self):
        with tempfile.TemporaryDirectory() as directory:
            fresh = Path(directory) / "fresh.json"
            fresh.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(fresh.read_bytes(), (audit.HERE / "evidence.json").read_bytes())


if __name__ == "__main__":
    unittest.main()
