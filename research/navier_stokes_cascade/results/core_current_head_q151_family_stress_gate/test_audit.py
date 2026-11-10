import json
import tempfile
import unittest
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import audit
import q151 as q


class Tests(unittest.TestCase):
    def q150(self):
        path = q.HERE.parent / "core_current_head_q150_initial_covariance_mean_join/evidence.json"
        return json.loads(path.read_text())["symbolic"]

    def test_conditional_canonical_feedback_is_closed(self):
        out = q.conditional_join(self.q150(), q.IDENTITY)
        self.assertIn("proved conditionally", out["canonical_feedback_inequality_D"])
        self.assertIn("axial covariance", out["canonical_postfast_scalars"])
        self.assertIn("A1", out["canonical_fresh_fast"])
        self.assertEqual(out["canonical_next_regular_source"], "D with pole 9p and additive tails excluded")

    def test_physical_family_is_not_transferred(self):
        out = q.conditional_join(self.q150(), q.IDENTITY)
        self.assertFalse(out["q150_inherited_physical_family_transferred_to_canonical_cycle"])
        self.assertIn("not the rebuilt canonical", out["reason"])

    def test_missing_provider_or_false_family_flag_rejected(self):
        original = q.provider_flags
        try:
            flags = deepcopy(original())
            flags["core_fast_row_gain_join"]["conditional_mean_feedback_D_and_next_regular_source_D_follow"] = False
            q.provider_flags = lambda: flags
            with self.assertRaises(ValueError):
                q.conditional_join(self.q150(), q.IDENTITY)
            flags = deepcopy(original())
            flags["canonical_primary_foundation"]["canonical_and_inherited_physical_primaries_or_roots_exactly_equal"] = True
            q.provider_flags = lambda: flags
            with self.assertRaises(ValueError):
                q.conditional_join(self.q150(), q.IDENTITY)
        finally:
            q.provider_flags = original

    def test_exact_rates(self):
        out = q.rates()
        self.assertEqual(out["B"], Q(99999, 100000))
        self.assertEqual(out["C"], Q(149999, 100000))
        self.assertEqual(out["D"], Q(74999, 50000))
        self.assertEqual(out["D_minus_1"], Q(24999, 50000))
        self.assertEqual(out["canonical_combined_feedback"], out["D"])

    def test_principal_stress_gate_and_closeness_obstruction(self):
        out = q.stress_gate(Q(1, 200))
        self.assertEqual(out["D_minus_1"], Q(24999, 50000))
        self.assertEqual(out["exponential_rate_in_ell"], Q(24999, 10000000))
        self.assertTrue(out["positive_rate"])
        self.assertFalse(out["polynomial_closeness_suffices"])
        self.assertIn("epsilon^-1 L_e", out["required_inequality"])

    def test_high_jet_request_is_finite_and_explicit(self):
        out = q.a1_request(0)
        self.assertEqual(out["compact_flatness_iterations"], 99999)
        self.assertEqual(out["old_prefast_source_order"], 200001)
        self.assertEqual(out["order_zero_fresh_fast_old_source"], 200006)
        self.assertFalse(out["actual_same_family_enclosure_certified_here"])

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError):
            q.conditional_join(self.q150(), "other")
        for k in (0, Q(1, 4), -1):
            with self.assertRaises(ValueError):
                q.rates(k)
        for h in (0, Q(1, 100), -1):
            with self.assertRaises(ValueError):
                q.stress_gate(h)
        with self.assertRaises(ValueError):
            q.a1_request(-1)

    def test_scope_record(self):
        record = q.record()
        self.assertFalse(record["inherited_and_canonical_primary_families_identified"])
        self.assertFalse(record["complete_residual_or_blowup_certified"])
        self.assertIn("O(1/L)", record["obstruction"])

    def test_frozen_record(self):
        with tempfile.TemporaryDirectory() as directory:
            fresh = Path(directory) / "fresh.json"
            fresh.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(fresh.read_bytes(), (audit.HERE / "evidence.json").read_bytes())


if __name__ == "__main__":
    unittest.main()
