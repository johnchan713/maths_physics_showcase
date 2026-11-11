import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q168


class Q168Tests(unittest.TestCase):
    def test_same_object_identity_chain(self):
        row = q168.provider_contracts()
        self.assertTrue(row["all_identity_links_hold"])
        self.assertTrue(all(row["identity_chain"].values()))
        self.assertTrue(row["same_physical_field_cutoff_sequence_support_and_pressure_gauge"])
        self.assertTrue(row["Q154_supplies_Q8_punctured_endpoint_premise"])
        self.assertTrue(row["Q155_compact_force_and_energy_conditional"])

    def test_exact_twelve_term_localization(self):
        row = q168.localization_contract()
        self.assertEqual(row["term_count"], 12)
        self.assertEqual(row["Q155_term_count"], 12)
        self.assertTrue(row["Q155_terms_match_after_notation_normalization"])
        self.assertEqual(row["complete_pressure"], "p0 grad c")
        self.assertEqual(row["unit_viscosity"], Q(1))
        self.assertTrue(row["force_divergence_retained"])
        for name in ("old_curl_transport", "curl_old_transport", "curl_square"):
            self.assertIn(name, row["term_names"])

    def test_only_thirteenth_term_vanishes_geometrically(self):
        row = q168.localization_contract()
        self.assertEqual(
            row["vanishing_thirteenth"],
            "(d dot grad c)u0=0 because d=grad(c) cross A",
        )
        self.assertNotIn("(d.grad c)u0", row["terms"])

    def test_zero_initial_divergence_free_localized_field(self):
        row = q168.localization_contract()
        self.assertTrue(row["divergence_free"])
        self.assertTrue(row["zero_initial_velocity"])
        self.assertTrue(row["protected_late_core_unchanged"])

    def test_derivative_request_and_reserve(self):
        row = q168.localization_contract()
        self.assertEqual(row["direct_force_orders"], {
            "potential": 7, "direct_swirl": 6, "pressure": 5,
        })
        self.assertEqual(row["Cauchy_orders"], {
            "potential": 8, "direct_swirl": 7, "pressure": 6,
        })
        self.assertEqual(row["Q167_supplied_orders"], {
            "potential": 9, "direct_swirl": 8, "pressure": 7,
        })
        self.assertEqual(set(row["unused_field_derivative"].values()), {1})
        self.assertEqual(row["residual_Cauchy_jet"], 6)

    def test_single_compact_force_diagonal(self):
        row = q168.force_contract()
        self.assertTrue(row["all_samples_have_geometric_tails"])
        self.assertFalse(row["raw_Taylor_convergence_used"])
        self.assertFalse(row["uniform_in_j_norm_assumed"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])
        for sample in row["sample_derivatives"]:
            self.assertGreater(sample["first_controlled_index"], sample["time_order"])

    def test_zero_safe_factor_two_energy_bound(self):
        row = q168.energy_contract(Q(7, 3))
        self.assertEqual(row["zero_safe_step"]["limit"], "sqrt(E(t))<=F(t)")
        self.assertEqual(row["sup_energy"], Q(49, 9))
        self.assertEqual(row["integrated_dissipation"], Q(49, 18))
        self.assertTrue(row["budgets_match"])
        self.assertTrue(row["factor_two_retained"])
        self.assertFalse(row["combined_endpoint_extrema_budget_claimed"])

    def test_energy_rejects_negative_force_integral(self):
        with self.assertRaises(ValueError):
            q168.energy_contract(-1)

    def test_frontier_moves_to_pressure_comparison(self):
        row = q168.remaining_gate()
        self.assertIn("closed-conditionally", row["U5_localization_force_energy"])
        self.assertIn("fixed-gauge pressure", row["next"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])
        self.assertFalse(row["comparison_or_core_exclusion_replayed_here"])
        self.assertEqual(q168.frontier()[-1], {
            "id": "U6", "status": "next-Q156-pressure-comparison-audit",
        })

    def test_scope_rejects_unconditional_promotion(self):
        row = q168.record()
        self.assertTrue(row["U5_localization_force_energy_closed_conditionally"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])
        self.assertFalse(row["comparison_core_exclusion_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_identity_rejected(self):
        for fn in (q168.provider_contracts, q168.localization_contract,
                   q168.force_contract, q168.remaining_gate, q168.frontier):
            with self.assertRaises(ValueError):
                fn("wrong")
        with self.assertRaises(ValueError):
            q168.energy_contract(identity="wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 407)
        self.assertGreater(record["closure"]["edges"], 1358)
        self.assertGreater(record["closure"]["source_bindings"], 5520)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
