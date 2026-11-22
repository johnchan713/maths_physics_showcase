import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q203


class Q203Tests(unittest.TestCase):
    def test_provider_statuses_and_same_object(self):
        self.assertTrue(all(q203.statuses().values()))
        self.assertTrue(all(q203.same_object().values()))

    def test_rate_route_correction(self):
        r = q203.rate_route()
        self.assertEqual(r["Q202_current_source_D"], Q(74999, 50000))
        self.assertEqual(r["Q202_cycle_gain"], Q(49997, 100000))
        self.assertEqual(r["log_margin"], Q(49997, 40000000))
        self.assertEqual(r["pure_margin"], Q(49996, 40000000))
        self.assertEqual(r["pure_margin_correction"], -Q(1, 40000000))
        self.assertFalse(r["Q202_D_may_be_read_as_same_exponent_pure_power"])
        self.assertTrue(r["pure_route_rebuild_and_scale_reselection_required"])

    def test_disjoint_compact_reserves(self):
        r = q203.rate_route()
        self.assertEqual(r["A0_own_margin"], Q(1, 40000000))
        self.assertEqual(r["A1_own_margin"], Q(1, 40000000))
        self.assertTrue(r["compact_accounts_disjoint"])

    def test_exact_diagonal_stage(self):
        log = q203.diagonal_stage()
        pure = q203.diagonal_stage(route="pure")
        self.assertEqual(log["J0"], 800)
        self.assertEqual(log["replay_stage"], 801)
        self.assertTrue(log["threshold_holds"])
        self.assertTrue(log["replay_beats_target_by_margin"])
        self.assertEqual(pure["J0"], 800)
        self.assertEqual(log["Q202_physical_common_pole"], 361200138)
        self.assertFalse(log["physical_pole_substituted_for_L_or_H"])

    def test_stage_with_nonzero_losses(self):
        r = q203.diagonal_stage(target=3, loss_L=5, loss_H=7)
        self.assertTrue(r["threshold_holds"])
        self.assertTrue(r["replay_beats_target_by_margin"])
        self.assertGreater(r["J0"], 0)

    def test_complete_six_group_diagonal(self):
        r = q203.complete_diagonal()
        self.assertEqual(len(r["six_groups"]), 6)
        self.assertTrue(r["all_six_groups_retain_complete_memory"])
        self.assertTrue(r["every_scale_cutoff_curl_and_pressure_derivative_retained"])
        self.assertTrue(r["divergence_free_term_by_term"])
        self.assertTrue(r["locally_finite"])
        self.assertFalse(r["raw_series_convergence_used"])
        self.assertFalse(r["stage_uniform_constant_or_J_of_q_used"])
        self.assertTrue(r["physical_domain_is_separate_from_diagonal_losses"])

    def test_corner_and_endpoint(self):
        r = q203.endpoint_join()
        self.assertEqual(r["jacobian_floor"], Q(99, 100))
        self.assertTrue(r["arbitrary_fixed_cartesian_power_closes"])
        self.assertFalse(r["residual_flatness_used_as_endpoint_field_bound"])
        self.assertTrue(r["punctured_endpoint_bound_supplied_separately"])
        self.assertTrue(r["full_exterior_pressure_convolution_retained"])
        self.assertGreater(r["pressure_S2_excess"], 0)
        self.assertFalse(r["exterior_time_Taylor_radius_positive"])
        self.assertFalse(r["terminal_origin_included"])

    def test_localization_force_and_energy(self):
        r = q203.localization_force()
        self.assertTrue(r["divergence_free"])
        self.assertTrue(r["zero_initial_velocity"])
        self.assertTrue(r["protected_inner_core_unchanged"])
        self.assertEqual(r["term_count"], 12)
        self.assertEqual(len(r["terms"]), 12)
        self.assertIn("=0", r["vanishing_thirteenth"])
        self.assertEqual(r["unit_viscosity"], 1)
        self.assertTrue(r["complete_pressure"])
        self.assertTrue(r["force_divergence_retained"])
        self.assertTrue(r["one_compact_smooth_force_conditionally"])
        self.assertFalse(r["raw_Taylor_convergence_used"])
        self.assertFalse(r["uniform_in_index_norm_assumed"])
        self.assertTrue(r["factor_two_retained"])
        self.assertFalse(r["future_velocity_or_pressure_constructed"])

    def test_precise_limit_and_scope(self):
        r = q203.obstruction()
        self.assertFalse(r["Q202_physical_pole_determines_L_H"])
        self.assertFalse(r["actual_L_H_or_giant_seminorms_evaluated"])
        self.assertLess(r["sourced_velocity_zero_jet_power"], 0)
        self.assertLess(r["sourced_pressure_zero_jet_power"], 0)
        self.assertLess(r["next_regular_source_zero_jet_power"], 0)
        self.assertFalse(r["available_upper_bounds_certify_absolute_decay"])
        review = q203.review()
        self.assertTrue(all(review["checks"].values()))
        self.assertFalse(review["scope"]["uniform_actual_pulse_or_complete_global_residual_bound_proved"])
        self.assertFalse(review["scope"]["unconditional_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q203.statuses, q203.same_object, q203.rate_route, q203.complete_diagonal,
                     q203.endpoint_join, q203.localization_force, q203.obstruction, q203.review):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q203.diagonal_stage(route="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 452)
        self.assertGreaterEqual(record["closure"]["edges"], 1631)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6093)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
