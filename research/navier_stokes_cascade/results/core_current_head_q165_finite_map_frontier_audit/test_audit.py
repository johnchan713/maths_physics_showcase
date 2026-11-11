import json
import tempfile
import unittest
from decimal import Decimal
from fractions import Fraction as Q
from pathlib import Path
import audit
import q165


class Q165Tests(unittest.TestCase):
    def test_Q2_regional_cover_and_corrected_ramp(self):
        row = q165.regional_contract()
        self.assertEqual(row["region_count"], 10)
        self.assertTrue(row["exact_shared_endpoints"])
        self.assertGreater(Decimal(row["corrected_outer_ramp_floor"]), Decimal("0.586"))
        self.assertTrue(row["floor_exceeds_one_half"])

    def test_Q2_flat_edge_obstruction_is_not_erased(self):
        row = q165.regional_contract()
        self.assertEqual(row["absolute_flat_edge_stress_infimum"], "0")
        self.assertFalse(row["positive_absolute_flat_edge_floor_claimed"])
        self.assertEqual(row["weighted_flat_edge_bound"], "|T0|>=c*zeta")
        self.assertTrue(row["smooth_edge_directions"])
        self.assertTrue(row["all_same_band_rows_positive"])
        self.assertFalse(row["uncertified_historical_existing_tilt_gate_used"])

    def test_Q3_frame_and_pressure_terms(self):
        row = q165.frame_contract()
        self.assertTrue(row["positive_normal_frame_and_diffusion_denominators"])
        self.assertTrue(row["picard_factorial_ratio_at_most_one_half"])
        self.assertTrue(row["unprojected_source_pressure_retained"])
        self.assertTrue(row["moving_normal_pressure_retained"])
        self.assertFalse(row["old_numeric_Cphi_gate_required"])
        self.assertFalse(row["unspecified_historical_first_band_certified"])
        self.assertTrue(row["potential_and_pressure_cutoff_commutators_retained"])

    def test_Q4_projectors_have_exact_moments_and_actual_supports(self):
        row = q165.projector_contract()
        self.assertGreater(row["fixed_order_majorant"], 0)
        self.assertEqual(len(row["projector_rows"]), 6)
        self.assertTrue(row["all_six_moments_exact"])
        self.assertTrue(row["actual_support_strictly_inside_hull"])
        self.assertTrue(row["crho_Z_times_P_retained"])
        self.assertTrue(row["same_choice_recompute_after_projector_selection"])
        self.assertFalse(row["uniformity_over_arbitrary_profiles_or_all_orders_claimed"])

    def test_Q5_log_reserve_corrects_zero_margin_pure_power_route(self):
        row = q165.weighted_inverse_contract()
        self.assertEqual(row["log_reserve"]["gain"], Q(1, 20000000))
        self.assertEqual(row["log_reserve"]["retained"], Q(1, 40000000))
        self.assertTrue(row["retained_log_reserve_positive"])
        self.assertTrue(row["old_pure_power_margin_is_exactly_zero"])
        self.assertTrue(row["physical_characteristic_shift_retained"])
        self.assertTrue(row["A1_axial_velocity_retained"])
        self.assertTrue(row["complete_moving_pressure_and_Pi1_Pi2_retained"])
        self.assertFalse(row["stage_or_derivative_uniform_constant_required"])

    def test_Q5_derivative_requests_remain_distinct_and_finite(self):
        requests = q165.weighted_inverse_contract()["request_graph"]
        self.assertEqual(requests["fast_tangential_or_axial"], 4)
        self.assertEqual(requests["fast_radial_potential"], 5)
        self.assertEqual(requests["whole_axial_pressure_moment"], 9)
        self.assertEqual(requests["A0_with_extra_integration"], 399997)
        self.assertEqual(requests["A1_with_extra_integration"], 200003)

    def test_current_head_provider_contracts_close_U2(self):
        row = q165.provider_contracts()
        self.assertTrue(row["Q164_U2_frontier_was_stale_relative_to_live_Q2_Q5_chain"])
        self.assertTrue(row["Q2_regional_and_factored_edge_contract_bound"])
        self.assertTrue(row["Q3_full_frame_and_whole_path_contract_bound"])
        self.assertTrue(row["Q4_fixed_projector_and_compensator_contract_bound"])
        self.assertTrue(row["Q5_weighted_inverse_and_log_source_contract_bound"])
        self.assertFalse(row["Q6_complete_residual_replayed_here"])

    def test_frontier_moves_to_Q6_residual(self):
        gate = q165.remaining_gate()
        self.assertIn("closed-by-current-head-Q2-Q5", gate["U2_finite_regional_frame_projector_inverse_source_chain"])
        self.assertIn("unit viscosity", gate["next"])
        self.assertFalse(gate["Q6_current_head_residual_audited_here"])
        self.assertFalse(gate["unconditional_theorem_gate_closed"])
        rows = q165.frontier()
        self.assertEqual(rows[-2], {"id": "U2", "status": "closed-by-Q165-finite-map-audit"})
        self.assertEqual(rows[-1], {"id": "U3", "status": "next-Q6-complete-residual-audit"})

    def test_scope_rejects_unconditional_promotion(self):
        row = q165.record()
        self.assertTrue(row["U2_finite_map_and_weighted_source_chain_closed"])
        self.assertFalse(row["Q6_complete_residual_audited_here"])
        self.assertFalse(row["unconditional_pulse_residual_force_exclusion_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_identity_rejected(self):
        for fn in (q165.regional_contract, q165.frame_contract,
                   q165.projector_contract, q165.weighted_inverse_contract,
                   q165.provider_contracts, q165.remaining_gate, q165.frontier):
            with self.assertRaises(ValueError):
                fn("wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 404)
        self.assertGreater(record["closure"]["edges"], 1343)
        self.assertGreater(record["closure"]["source_bindings"], 5495)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
