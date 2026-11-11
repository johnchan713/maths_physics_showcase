import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q166


class Q166Tests(unittest.TestCase):
    def test_Q165_and_Q6_share_the_identical_Q5_object(self):
        row = q166.provider_contracts()
        self.assertTrue(row["Q165_U2_closed"])
        self.assertTrue(row["Q165_frontier_was_Q6"])
        self.assertTrue(row["identical_Q5_evidence_enters_Q165_and_Q6"])
        self.assertTrue(row["Q6_all_provider_gates_hold"])
        self.assertFalse(row["Q6_unconditional_claimed"])

    def test_universal_selected_parameter_margins(self):
        row = q166.universal_margin()
        self.assertEqual(row["delta"], Q(49997, 100000))
        self.assertEqual(row["one_stage_gain"], Q(49997, 20000000))
        self.assertEqual(row["regular_retained_margin"], Q(49997, 40000000))
        self.assertEqual(row["limiting_cross_margin_at_least"], Q(49997, 40000000))
        self.assertEqual(row["A0_own_log_retained_margin"], Q(1, 40000000))
        self.assertEqual(row["A1_own_log_retained_margin"], Q(1, 40000000))
        self.assertTrue(row["regular_and_cross_reserves_same_size_but_distinct_owners"])
        self.assertFalse(row["compact_reserves_may_pay_regular_or_cutoff_logs"])
        self.assertFalse(row["disjoint_reserves_added_to_strengthen_common_exponent"])

    def test_log_absorption_requires_positive_owned_reserve(self):
        self.assertEqual(q166.logarithmic_absorption(37, Q(1, 40000000)), 37**2 * 10000000)
        with self.assertRaises(ValueError):
            q166.logarithmic_absorption(-1, Q(1, 2))
        with self.assertRaises(ValueError):
            q166.logarithmic_absorption(1, 0)

    def test_complete_six_group_samples(self):
        row = q166.residual_contract()
        self.assertEqual(row["sample_stages"], [12817, 18418, 37619])
        self.assertTrue(row["all_sample_controls_hold"])
        self.assertTrue(row["all_six_groups_retained_in_order"])
        self.assertTrue(row["both_cross_orientations_retained"])
        self.assertTrue(row["complete_pressure_covariance_axial_and_cutoff_memory_retained"])

    def test_default_minimum_is_the_cross_not_the_fixed_state(self):
        row = q166.residual_contract()
        self.assertEqual(Q(row["default_limiting_cross_exponent"]), Q(20030773, 20000000))
        self.assertEqual(row["default_complete_minimum"], row["default_limiting_cross_exponent"])
        self.assertLess(Q(row["default_complete_minimum"]), Q(row["default_fixed_state_exponent"]))
        self.assertTrue(row["complete_exponent_is_minimum_not_sum"])

    def test_fixed_jet_then_finite_stage_quantifier(self):
        row = q166.residual_contract()
        self.assertTrue(row["fixed_jet_target_then_finite_stage_quantifier_retained"])
        self.assertFalse(row["stage_uniform_constant_or_J_of_q_used"])

    def test_frontier_moves_to_Q7_without_downstream_promotion(self):
        gate = q166.remaining_gate()
        self.assertIn("closed-conditionally", gate["U3_complete_six_group_local_residual"])
        self.assertIn("Q7", gate["next"])
        self.assertFalse(gate["closed_corner_endpoint_localization_force_energy_comparison_replayed_here"])
        self.assertFalse(gate["unconditional_theorem_gate_closed"])
        self.assertEqual(q166.frontier()[-1], {"id": "U4", "status": "next-Q7-corner-and-endpoint-audit"})

    def test_scope_rejects_unconditional_promotion(self):
        row = q166.record()
        self.assertTrue(row["U3_complete_local_residual_closed_conditionally"])
        self.assertFalse(row["corner_endpoint_or_localized_force_replayed_here"])
        self.assertFalse(row["unconditional_pulse_residual_force_exclusion_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_identity_rejected(self):
        for fn in (q166.provider_contracts, q166.universal_margin,
                   q166.residual_contract, q166.remaining_gate, q166.frontier):
            with self.assertRaises(ValueError):
                fn("wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 405)
        self.assertGreater(record["closure"]["edges"], 1349)
        self.assertGreater(record["closure"]["source_bindings"], 5505)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
