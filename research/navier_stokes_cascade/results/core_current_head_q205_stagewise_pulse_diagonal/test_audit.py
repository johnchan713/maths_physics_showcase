import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q205


class Q205Tests(unittest.TestCase):
    def test_provider_statuses_and_same_object(self):
        self.assertTrue(all(q205.statuses().values()))
        self.assertTrue(all(q205.same_object_links().values()))

    def test_fixed_cycle_actual_bounds(self):
        r = q205.fixed_cycle_bounds()
        self.assertEqual(r["velocity_zero_jet_power"], Q(-4975001, 10000000))
        self.assertEqual(r["pressure_zero_jet_power"], Q(-10000001, 10000000))
        self.assertTrue(r["constant_independent_of_running_band_at_fixed_cycle"])
        self.assertTrue(r["constant_may_depend_on_cycle_derivative_and_profiles"])
        self.assertTrue(r["full_path_zero_gives_smooth_zero_extension"])
        self.assertTrue(r["compact_source_need_not_have_zero_exit"])
        self.assertTrue(r["potential_and_pressure_cut_before_curl"])
        self.assertTrue(r["closed_support_not_enlarged_by_curl_graph_or_frame"])
        self.assertTrue(r["moving_normal_pressure_retained"])
        self.assertTrue(r["source_normal_pressure_retained"])

    def test_uniform_raw_countermodel(self):
        r = q205.raw_uniform_obstruction(count=20)
        self.assertEqual(r["countermodel_Cj"], "2^(j^2)")
        self.assertTrue(r["countermodel_satisfies_stagewise_representative_bound"])
        self.assertEqual(r["fixed_q0"], Q(1, 2))
        self.assertTrue(r["strictly_increasing_sample"])
        self.assertTrue(r["unbounded_in_j"])
        self.assertFalse(r["common_raw_hierarchy_constant_follows"])
        self.assertFalse(r["actual_pulse_lower_bound"])
        self.assertFalse(r["residual_lower_bound"])

    def test_diagonal_absorbs_arbitrary_stage_constant(self):
        for j, p in ((1, 0), (10, 7), (50, 19)):
            r = q205.diagonal_absorption(j=j, P=p)
            self.assertTrue(r["endpoint_bound_holds"])
            self.assertTrue(r["log_bound_decreasing_after_cutoff"])
            self.assertTrue(r["works_for_arbitrary_finite_stage_constants"])
            self.assertTrue(r["scale_chosen_after_current_constants"])
            self.assertFalse(r["stage_uniform_constant_used"])
            self.assertFalse(r["J_of_q_substitution_used"])
            self.assertTrue(r["locally_finite_on_q_positive"])

    def test_pressure_support_memory(self):
        r = q205.pressure_support_memory()
        self.assertTrue(r["cut_complete_potential_direct_swirl_and_pressure"])
        self.assertFalse(r["stage_dependent_potential_gauge_allowed"])
        self.assertFalse(r["stage_dependent_pressure_gauge_allowed"])
        self.assertTrue(r["six_covariances"])
        self.assertTrue(r["radial_and_axial_velocity"])
        self.assertTrue(r["moving_pressure_product"])
        self.assertTrue(r["Pi1_Pi2_and_same_base_tails"])
        self.assertTrue(r["viscosity_cutoffs_both_crosses_and_square"])
        self.assertTrue(r["reserved_pressure_floor_preserved"])

    def test_review_scope(self):
        r = q205.review()
        self.assertTrue(all(r["checks"].values()))
        self.assertTrue(r["scope"]["stagewise_actual_generated_pulse_bounds_joined"])
        self.assertFalse(r["scope"]["uniform_generated_hierarchy_constant_required_by_diagonal"])
        self.assertFalse(r["scope"]["uniform_uncut_generated_hierarchy_proved"])
        self.assertFalse(r["scope"]["absolute_pulse_decay_proved"])
        self.assertFalse(r["scope"]["unconditional_complete_residual_or_blowup"])

    def test_invalid_identity_and_arguments(self):
        for call in (q205.statuses, q205.same_object_links, q205.fixed_cycle_bounds,
                     q205.raw_uniform_obstruction, q205.diagonal_absorption,
                     q205.pressure_support_memory, q205.review, q205.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q205.raw_uniform_obstruction(count=1)
        with self.assertRaises(ValueError):
            q205.diagonal_absorption(j=0)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 454)
        self.assertGreaterEqual(record["closure"]["edges"], 1644)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6114)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
