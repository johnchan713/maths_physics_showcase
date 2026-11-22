import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q204


class Q204Tests(unittest.TestCase):
    def test_provider_statuses_and_same_object_chain(self):
        self.assertTrue(all(q204.statuses().values()))
        self.assertTrue(all(q204.same_object_chain().values()))

    def test_global_pressure_moment_ledger(self):
        r = q204.global_pressure_moments()
        self.assertFalse(r["raw_I_infinity_exists"])
        self.assertEqual(r["global_zero_defects"], ("M", "R", "J", "S", "Pi_axis+Cp"))
        self.assertTrue(r["functional_in_eta"])
        self.assertTrue(r["intermediate_pressure_memory_may_be_nonzero"])
        self.assertTrue(r["same_pressure_gauge_gives_backward_tail"])
        self.assertTrue(r["zero_five_defects_remove_all_exterior_memory"])
        self.assertTrue(r["pressure_eta_slope_defect_grows_linearly_in_X"])
        self.assertFalse(r["zero_U_removes_M_induced_radial_velocity"])
        self.assertTrue(r["all_exact_checks"])

    def test_coefficient_supports(self):
        r = q204.coefficient_supports()
        self.assertEqual(r["inverse_bound"], 800)
        self.assertEqual(r["order_1_support"], {"velocity": "Xvel", "pressure": "Xv", "stress": "Xb"})
        self.assertEqual(r["order_2_support"], {"velocity": "Xvel", "pressure": "Xv", "stress": "Xv"})
        self.assertEqual(set(r["order_3_support"].values()), {"Xvel"})
        self.assertEqual(r["reserved_pressure_floor_sample"], Q(9, 1400))
        self.assertTrue(r["pressure_floor_independent_of_K1_K2"])
        self.assertTrue(r["pressure_floor_sharp_after_constant_cancellation"])

    def test_rule_schedule_and_uniform_primary(self):
        r = q204.schedule_and_primary()
        self.assertEqual(r["schedule_binary_margin_sample"], 2)
        self.assertTrue(r["schedule_strict_padded_rule_holds"])
        self.assertTrue(r["scale_chosen_after_finite_Cnm"])
        self.assertFalse(r["preexisting_numeric_scale_certificate_required"])
        self.assertTrue(r["uniform_primary_pulse_family_certified"])
        self.assertFalse(r["canonical_and_physical_families_identified"])
        self.assertFalse(r["complete_generated_hierarchy_certified"])
        self.assertEqual(r["primary_zero_jet_powers"]["velocity"], Q(-201, 400))
        self.assertEqual(r["primary_zero_jet_powers"]["pressure"], Q(-201, 200))

    def test_complete_covariance_and_current_cycle(self):
        r = q204.covariance_and_current_cycle()
        self.assertTrue(r["same_physical_family"])
        self.assertTrue(r["principal_stress_gate_closed_exactly"])
        self.assertEqual(r["complete_curl_margin"], Q(1, 2))
        self.assertEqual(r["current_source_order"], Q(74999, 50000))
        self.assertTrue(r["current_source_has_complete_operator_memory"])
        self.assertTrue(r["moving_and_source_normal_pressure_retained"])
        self.assertEqual(r["next_source"], Q(39999, 20000))
        self.assertEqual(r["next_target"], Q(49999, 20000))
        self.assertEqual(r["target_slack"], Q(1, 50000))
        self.assertTrue(r["cycle_closes_relative_to_named_inputs"])

    def test_complete_diagonal_and_localized_residual(self):
        r = q204.diagonal_residual()
        self.assertEqual(len(r["all_six_groups"]), 6)
        self.assertTrue(r["complete_memory"])
        self.assertTrue(r["cutoff_curl_and_pressure_derivatives"])
        self.assertTrue(r["fixed_jet_then_finite_stage"])
        self.assertFalse(r["stage_uniform_constant_used"])
        self.assertTrue(r["conditional_fixed_jet_flatness"])
        self.assertEqual(r["localization_term_count"], 12)
        self.assertTrue(r["complete_localized_pressure"])
        self.assertEqual(r["unit_viscosity"], 1)

    def test_precise_limits(self):
        r = q204.precise_limits()
        self.assertTrue(r["reserved_pressure_cannot_be_replaced_by_leading_pressure"])
        self.assertLess(r["primary_velocity_zero_jet_power"], 0)
        self.assertLess(r["primary_pressure_zero_jet_power"], 0)
        self.assertEqual(r["sourced_velocity_zero_jet_power"], Q(-4975001, 10000000))
        self.assertEqual(r["sourced_pressure_zero_jet_power"], Q(-10000001, 10000000))
        self.assertFalse(r["absolute_decay_certified"])
        self.assertTrue(r["uniform_primary_scale_bounds_certified"])
        self.assertFalse(r["uniform_generated_hierarchy_constants_certified"])
        self.assertFalse(r["Q202_constants_uniform_in_cycle"])

    def test_review_scope(self):
        r = q204.review()
        self.assertTrue(all(r["checks"].values()))
        self.assertTrue(r["scope"]["global_ledgers_and_uniform_primary_joined_to_Q203_conditionally"])
        self.assertFalse(r["scope"]["uniform_complete_generated_hierarchy_proved"])
        self.assertFalse(r["scope"]["unconditional_complete_residual_force_exclusion_or_blowup"])

    def test_invalid_identity(self):
        for call in (q204.statuses, q204.same_object_chain, q204.global_pressure_moments,
                     q204.coefficient_supports, q204.schedule_and_primary,
                     q204.covariance_and_current_cycle, q204.diagonal_residual,
                     q204.precise_limits, q204.review):
            with self.assertRaises(ValueError):
                call(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 453)
        self.assertGreaterEqual(record["closure"]["edges"], 1641)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6107)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
