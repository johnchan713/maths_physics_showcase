import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q163


class Q163Tests(unittest.TestCase):
    def test_pressure_recovery_and_force_cancellation(self):
        row = q163.comparison_contract()
        self.assertTrue(row["common_force_cancels_once"])
        self.assertTrue(row["finite_energy_essential"])
        self.assertFalse(row["competitor_pressure_decay_assumed"])
        self.assertIn("R_i R_j", row["pressure_recovery"])
        self.assertIn("-Delta pi", row["pressure_difference"])

    def test_expanding_ball_remainder_is_radius_minus_two(self):
        row = q163.comparison_contract()
        self.assertEqual(row["worst_radius_decay"], Q(2))
        self.assertEqual(row["compact_mixed_flux"], "O(R^-5/2)")
        self.assertEqual(row["cutoff_power"], "chi_R=phi_R^8")
        self.assertIn("R^-2", row["comparison"]["localized"])

    def test_comparison_has_only_ordinary_competitor_hypotheses(self):
        row = q163.comparison_contract()["comparison"]
        self.assertFalse(row["competitor_derivative_bound_assumed"])
        self.assertFalse(row["competitor_boundary_flux_bound_assumed"])
        self.assertFalse(row["competitor_dissipation_bound_assumed"])
        self.assertFalse(row["terminal_uniform_constants_assumed"])

    def test_same_physical_core_parameters_and_support_hole(self):
        row = q163.core_contract()
        self.assertEqual(row["parameters"]["A"], Q(101, 200))
        self.assertEqual(row["parameters"]["X_in"], Q(1, 101))
        self.assertEqual(row["parameters"]["X_a"], Q(4, 101))
        self.assertTrue(row["path_inside_annular_support_hole"])
        self.assertTrue(row["localization_plateau_preserved"])
        self.assertFalse(row["binding"]["canonical_physical_primary_identification_used"])

    def test_combined_core_threshold_is_strictly_positive(self):
        row = q163.core_contract()
        self.assertTrue(row["threshold_positive"])
        self.assertEqual(row["threshold_example"]["q_core"], Q(1, 4) ** 100)
        self.assertEqual(row["threshold_example"]["plateau"], Q(1, 32))
        self.assertEqual(row["threshold_example"]["tau_star"], Q(1, 4) ** 100)

    def test_limit_order_uses_T_tau_between_path_time_and_one(self):
        row = q163.limit_order()
        self.assertTrue(row["strict_time_order"])
        self.assertEqual(row["T_tau"], 1 - row["tau"] / 2)
        self.assertEqual(row["t"], 1 - row["tau"])
        self.assertFalse(row["uniform_comparison_constant_as_T_to_one_required"])
        self.assertIn("R to infinity", row["second"])
        self.assertIn("tau", row["third"])

    def test_bad_protected_tau_rejected(self):
        threshold = q163.core_contract()["threshold_example"]["tau_star"]
        for tau in (0, -1, threshold * 2):
            with self.assertRaises(ValueError):
                q163.limit_order(tau)

    def test_conditional_contradiction_scope(self):
        row = q163.contradiction_contract()
        self.assertTrue(row["comparison_applies_at_path_time"])
        self.assertTrue(row["positive_exponent"])
        self.assertIn("no globally smooth", row["conditional_conclusion"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])
        self.assertFalse(row["unconditional_blowup_claimed"])

    def test_provider_contracts_close_P7_conditionally(self):
        row = q163.provider_contracts()
        self.assertTrue(row["same_Q162_field_force_pressure_gauge_and_localization"])
        self.assertTrue(row["fixed_T_then_radius_then_core_limit_order"])
        self.assertTrue(row["P7_conditional_consequence_closed"])
        self.assertTrue(row["historical_conditional_records_preserved"])

    def test_remaining_foundation_gate_is_not_promoted(self):
        row = q163.remaining_foundation_gate()
        self.assertFalse(row["comparison_foundation_all_dependencies_discharged"])
        self.assertFalse(row["Q1_every_provider_theorem_independently_reproved"])
        self.assertIn("lambda<1", row["required_inequalities"])
        self.assertIn("Pi_axis+C_p,infinity=0", row["required_inequalities"])
        self.assertFalse(row["dependency_hashes_or_passing_tests_suffice"])
        self.assertFalse(row["unconditional_theorem_gate_closed"])

    def test_frontier_moves_to_independent_foundation_audit(self):
        rows = q163.frontier()
        self.assertEqual([row["id"] for row in rows[:7]], [f"P{i}" for i in range(1, 8)])
        self.assertEqual(rows[6]["status"], "closed-as-conditional-consequence-by-Q163")
        self.assertEqual(rows[7], {"id": "U1", "status": "next-independent-foundation-audit"})

    def test_scope_rejects_unconditional_promotion(self):
        row = q163.record()
        self.assertTrue(row["P7_comparison_and_core_exclusion_replayed_conditionally"])
        self.assertFalse(row["unconditional_admissible_force_breakdown_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_identity_rejected(self):
        for fn in (q163.comparison_contract, q163.core_contract,
                   q163.limit_order, q163.contradiction_contract,
                   q163.provider_contracts, q163.remaining_foundation_gate,
                   q163.frontier):
            with self.assertRaises(ValueError):
                if fn is q163.limit_order:
                    fn(identity="wrong")
                else:
                    fn("wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 402)
        self.assertGreater(record["closure"]["edges"], 1324)
        self.assertGreater(record["closure"]["source_bindings"], 5468)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
