import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q170


class Q170Tests(unittest.TestCase):
    def test_same_object_identity_chain(self):
        row = q170.provider_contracts()
        self.assertTrue(row["all_identity_links_hold"])
        self.assertTrue(all(row["identity_chain"].values()))
        self.assertTrue(row["same_finalized_background"])
        self.assertTrue(row["physical_diagonal_has_unchanged_annular_support_hole"])
        self.assertTrue(row["localization_plateau_bound_to_same_field"])
        self.assertTrue(row["Q169_preterminal_comparison"])

    def test_exact_path_and_support_hole(self):
        row = q170.core_path_contract(Q(1, 1000))
        self.assertEqual(row["h"], Q(1, 200))
        self.assertEqual(row["A"], Q(101, 200))
        self.assertEqual(row["X_in"], Q(1, 101))
        self.assertEqual(row["X_a"], Q(4, 101))
        self.assertEqual(row["support_hole_margin"], Q(3, 101))
        self.assertEqual(row["coordinates"]["q"], Q(1, 1000))
        self.assertEqual(row["coordinates"]["eta"], 0)
        self.assertEqual(row["coordinates"]["X"], Q(1, 101))

    def test_background_smallness_is_not_annular_vanishing(self):
        row = q170.core_path_contract()
        self.assertFalse(row["background_terms_vanish_exactly"])
        self.assertTrue(row["annular_diagonal_terms_vanish_exactly"])
        self.assertIn("q^(2h)", row["background_remainder"])
        self.assertFalse(row["canonical_physical_primary_identification_used"])
        self.assertFalse(row["pressure_gauge_changed"])

    def test_exponent_one_hundred_threshold(self):
        row = q170.threshold_contract()
        self.assertEqual(row["q_core"], Q(1, 4) ** 100)
        self.assertEqual(row["plateau"], Q(1, 32))
        self.assertEqual(row["tau_star"], Q(1, 4) ** 100)
        self.assertTrue(row["strictly_positive"])
        self.assertTrue(row["normalized_error_within_half"])
        self.assertEqual(row["normalized_error_bound_on_interval"], Q(1, 2))

    def test_localization_plateau_is_quantitative(self):
        row = q170.threshold_contract()
        self.assertTrue(row["time_plateau_bound"])
        self.assertTrue(row["radial_plateau_bound"])
        self.assertTrue(row["path_inside_quarter_radius_squared"])
        self.assertIn("101/200", row["lower_bound"])

    def test_fixed_T_radius_then_terminal_order(self):
        row = q170.limit_contract()
        self.assertTrue(row["strict_time_order"])
        self.assertEqual(row["path_time"], 1 - row["tau"])
        self.assertEqual(row["comparison_time"], 1 - row["tau"] / 2)
        self.assertIn("R to infinity", row["second"])
        self.assertIn("tau", row["third"])
        self.assertFalse(row["uniform_comparison_constant_as_T_to_one_required"])

    def test_conditional_competitor_exclusion_scope(self):
        row = q170.limit_contract()
        self.assertEqual(row["growth_exponent"], Q(101, 200))
        self.assertIn("no globally smooth", row["conditional_conclusion"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])
        self.assertFalse(row["unconditional_blowup_claimed"])

    def test_bad_tau_is_rejected(self):
        threshold = q170.threshold_contract()["tau_star"]
        for tau in (0, -1, threshold * 2):
            with self.assertRaises(ValueError):
                q170.limit_contract(tau)

    def test_frontier_moves_to_provider_discharge_audit(self):
        row = q170.remaining_gate()
        self.assertIn("closed-conditionally", row["U7_protected_core_and_limit_order"])
        self.assertTrue(row["Q163_stale_next_core_audit_already_completed_by_Q164"])
        self.assertFalse(row["all_named_continuum_providers_independently_discharged_here"])
        self.assertFalse(row["unconditional_admissible_force_or_blowup_certified"])
        self.assertEqual(q170.frontier()[-1], {
            "id": "U8", "status": "next-consolidated-provider-discharge-audit",
        })

    def test_scope_rejects_unconditional_promotion(self):
        row = q170.record()
        self.assertTrue(row["U7_protected_core_limit_closed_conditionally"])
        self.assertTrue(row["conditional_same_force_competitor_exclusion"])
        self.assertFalse(row["unconditional_admissible_force_or_blowup_certified"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_identity_rejected(self):
        for fn in (q170.provider_contracts, q170.core_path_contract,
                   q170.threshold_contract, q170.remaining_gate,
                   q170.frontier):
            with self.assertRaises(ValueError):
                fn(identity="wrong")
        with self.assertRaises(ValueError):
            q170.limit_contract(identity="wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 409)
        self.assertGreater(record["closure"]["edges"], 1368)
        self.assertGreater(record["closure"]["source_bindings"], 5540)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
