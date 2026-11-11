import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q164


class Q164Tests(unittest.TestCase):
    def test_all_index_operator_contract(self):
        row = q164.operator_contract()
        self.assertEqual(row["product_constant_A0"], 256)
        self.assertEqual(row["integrated_derivative_constant_D0"], "80/rho")
        self.assertTrue(row["all_index_product_bound_rederived"])
        self.assertTrue(row["all_index_mixed_allocation_bound_rederived"])
        self.assertTrue(all(row["exact_checks"].values()))
        self.assertTrue(row["finite_samples_are_not_all_index_proofs"])

    def test_complete_21_term_pressure_ledger(self):
        row = q164.source_contract()
        self.assertEqual(row["term_count"], 21)
        self.assertEqual(row["unique_term_names"], 21)
        self.assertEqual(row["group_sums"]["R1_linear"], Q(19237, 50))
        self.assertEqual(row["group_sums"]["R1_quadratic"], Q(952, 25))
        self.assertEqual(row["group_sums"]["R2_pressure"], 34)
        self.assertEqual(row["pressure_primitive"], 80 * 256**3)
        self.assertEqual(row["pressure_source"], 2720 * 256**4)
        self.assertEqual(row["pressure_lipschitz_per_S"], 5440 * 256**4)
        self.assertTrue(row["complete_pressure_memory_retained"])

    def test_selected_parameter_formula_is_strict(self):
        row = q164.contraction_contract()
        self.assertEqual(row["epsilon"], Q(1, 10**16))
        self.assertEqual(row["coefficient_before_H"], Q(1, 10**50))
        self.assertEqual(row["H_denominator_power"], 127)
        self.assertTrue(row["strictly_below_ceiling_because_H_gt_one"])
        self.assertTrue(row["ball_invariant"])
        self.assertTrue(row["strict_contraction"])
        self.assertFalse(row["finite_jet_iteration_used_as_convergence_proof"])

    def test_original_equations_and_regular_forward_stress(self):
        row = q164.equation_contract()
        self.assertIn("Pi_X=phi^2/C^2", row["pressure_identity"])
        self.assertTrue(row["axis_regularity_forces_both_constants_zero"])
        self.assertEqual(row["exact_leading_core_stress"], "T0 identically zero")
        self.assertFalse(row["pressure_or_one_of_its_derivatives_removed"])

    def test_same_gauge_leading_pressure_only(self):
        row = q164.pressure_contract()
        self.assertEqual(row["A"], Q(101, 200))
        self.assertTrue(row["positive_exponent"])
        self.assertFalse(row["nonzero_mode_endpoint_bounded"])
        self.assertTrue(row["smooth_endpoint_compatibility_forces_C_zero"])
        self.assertTrue(row["leading_core_equals_completed_leading_backward_pressure"])
        self.assertFalse(row["full_corrected_pressure_collapsed_to_leading_pressure"])
        self.assertTrue(row["later_background_wave_mean_localization_pressure_retained"])

    def test_current_head_provider_contracts_correct_q163_frontier(self):
        row = q164.provider_contracts()
        self.assertTrue(row["Q163_U1_was_stale_relative_to_live_provider_graph"])
        self.assertTrue(row["independent_all_index_rederivation_bound_to_Q1_current_head"])
        self.assertTrue(row["original_equations_and_regular_forward_stress_recovered"])
        self.assertTrue(row["same_gauge_leading_pressure_identified"])
        self.assertTrue(row["later_complete_pressure_memory_retained"])

    def test_frontier_moves_to_finite_maps(self):
        gate = q164.remaining_gate()
        self.assertIn("closed-by-existing-independent", gate["U1_all_index_core_contraction_original_equations_leading_pressure"])
        self.assertIn("regional cover", gate["next"])
        self.assertFalse(gate["Q2_regional_continuum_infima_discharged"])
        self.assertFalse(gate["Q3_full_frame_primitive_domain_discharged"])
        self.assertFalse(gate["dependency_hashes_or_passing_tests_suffice"])
        rows = q164.frontier()
        self.assertEqual(rows[-2], {"id": "U1", "status": "closed-by-independent-core-pressure-audit"})
        self.assertEqual(rows[-1], {"id": "U2", "status": "next-independent-finite-map-audit"})

    def test_scope_rejects_unconditional_promotion(self):
        row = q164.record()
        self.assertTrue(row["U1_core_contraction_and_leading_pressure_closed"])
        self.assertFalse(row["Q2_or_Q3_independently_closed"])
        self.assertFalse(row["unconditional_pulse_residual_force_exclusion_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_identity_rejected(self):
        for fn in (q164.operator_contract, q164.source_contract,
                   q164.contraction_contract, q164.equation_contract,
                   q164.pressure_contract, q164.provider_contracts,
                   q164.remaining_gate, q164.frontier):
            with self.assertRaises(ValueError):
                fn("wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 403)
        self.assertGreater(record["closure"]["edges"], 1337)
        self.assertGreater(record["closure"]["source_bindings"], 5485)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
