import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q162


class Q162Tests(unittest.TestCase):
    def test_complete_diagonal_contract(self):
        row = q162.diagonal_contract()
        self.assertTrue(row["all_fixed_graphs_finite_and_acyclic"])
        self.assertTrue(row["complete_pressure_covariance_axial_velocity_and_supports_retained"])
        self.assertTrue(row["arbitrary_finite_stage_constants_absorbed"])
        self.assertTrue(row["cutoff_sequence_chosen_once"])
        self.assertTrue(row["comparison_stage_fixed_before_q_limit"])
        self.assertTrue(row["one_locally_finite_divergence_free_field"])
        self.assertTrue(row["all_fixed_residual_jets_flat"])
        self.assertFalse(row["raw_series_convergence_used"])

    def test_complete_residual_groups(self):
        groups = q162.diagonal_contract()["complete_residual_groups"]
        for term in ("time", "viscosity", "pressure", "old_new", "new_old", "new_square"):
            self.assertIn(term, groups)

    def test_punctured_endpoint_window_and_margin(self):
        row = q162.endpoint_contract()
        self.assertEqual(row["implicit_floor"], Q(99, 100))
        self.assertEqual(row["closed_band_window"]["bands"], tuple(range(18, 27)))
        self.assertTrue(row["observation_inside_qstar"])
        self.assertTrue(row["observation_beyond_first_band"])
        self.assertFalse(row["singular_origin_included"])

    def test_pressure_memory_is_full_convolution(self):
        row = q162.endpoint_contract()
        self.assertTrue(row["pressure_convolution_strictly_exceeds_two_A2"])
        self.assertGreater(row["pressure_convolution_at_two"], 0)
        self.assertFalse(row["exterior_time_Taylor_radius_positive"])

    def test_endpoint_derivative_ledger(self):
        req = q162.endpoint_contract()["fixed_request_example"]
        self.assertEqual(req, {"m": 5, "potential": 9, "direct_swirl": 8,
                               "pressure": 7, "corner_residual": 6})

    def test_exact_twelve_term_localization(self):
        row = q162.localization_force_contract()
        self.assertEqual(row["term_count"], 12)
        self.assertEqual(len(row["localization_terms"]), 12)
        self.assertEqual(row["unit_viscosity"], 1)
        self.assertEqual(row["complete_pressure"], "p0 grad c")
        self.assertTrue(row["force_divergence_retained"])
        self.assertIn("=0", row["vanishing_thirteenth"])

    def test_force_derivative_and_corner_ledgers(self):
        row = q162.localization_force_contract()
        request = row["derivative_request_example"]
        self.assertEqual(request["force"], {"potential": 7, "direct_swirl": 6,
                                            "pressure": 5})
        self.assertEqual(request["Q154_spare"], {"potential": 1,
                                                 "direct_swirl": 1,
                                                 "pressure": 1})
        self.assertTrue(row["corner_transfer_arbitrary_power"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])

    def test_energy_arithmetic_and_scope(self):
        energy = q162.localization_force_contract()["energy"]
        self.assertEqual(energy["sup_energy"], Q(49, 9))
        self.assertEqual(energy["integrated_dissipation"], Q(49, 18))
        self.assertFalse(energy["sum_of_separate_extrema_budget_claimed"])

    def test_same_object_provider_contract(self):
        row = q162.provider_contracts()
        self.assertTrue(row["same_physical_field_cutoff_sequence_support_and_pressure_gauge"])
        self.assertTrue(row["P6_complete_diagonal_endpoint_localized_force_contract_closed"])
        self.assertTrue(row["historical_conditional_records_preserved"])
        self.assertTrue(row["P7_comparison_and_core_exclusion_not_replayed_here"])

    def test_frontier_advances_only_to_P7(self):
        rows = q162.frontier()
        self.assertEqual([row["id"] for row in rows], [f"P{i}" for i in range(1, 8)])
        self.assertEqual(rows[5]["status"], "closed-by-Q162-same-object-replay")
        self.assertEqual(rows[6]["status"], "next-replay-required")

    def test_scope_rejects_promotion(self):
        row = q162.record()
        self.assertTrue(row["P6_complete_diagonal_endpoint_localized_force_contract_discharged"])
        self.assertFalse(row["P7_comparison_and_core_exclusion_replayed"])
        self.assertFalse(row["unconditional_competitor_exclusion_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_identity_rejected(self):
        for fn in (q162.diagonal_contract, q162.endpoint_contract,
                   q162.localization_force_contract, q162.provider_contracts,
                   q162.frontier):
            with self.assertRaises(ValueError):
                fn("wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 400)
        self.assertGreater(record["closure"]["edges"], 1310)
        self.assertGreater(record["closure"]["source_bindings"], 5450)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
