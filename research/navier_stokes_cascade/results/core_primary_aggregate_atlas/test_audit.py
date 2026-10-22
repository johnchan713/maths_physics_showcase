#!/usr/bin/env python3
"""Actual atlas aggregate, edge constant and reference-convention tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import P,aggregate,atlas,edge_constant,rescale,review


class PrimaryAggregateAtlasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_actual_closed_support_separation_is_replayed(self):
        self.controls("actual_atlas","closed_enlarged_absolute_preimages_are_disjoint")
        self.controls("actual_atlas","same_band_closed_doubled_boxes")
        self.assertEqual(atlas.radius_bound(),Q(1,6057072))

    def test_actual_covering_gap_and_raw_chain_factor(self):
        self.controls("actual_atlas","all_band_half_line_gap_bound_below_three")
        self.controls("actual_atlas","symmetric_covering_has_same_column_norm")
        self.controls("actual_atlas","raw_jet_pullback_bound")

    def test_deck_average_is_normalized_without_a_copy_factor(self):
        self.controls("actual_atlas","normalized_deck_average_has_no_degree_factor")
        self.controls("actual_atlas","unnormalized_copy_sum_is_different")
        self.assertEqual(atlas.counting_review()["common_lift_factor"],2744)

    def test_slow_label_count_is_separate_from_actual_auxiliary_overlap(self):
        self.assertEqual(atlas.counting_review()["original_pointwise_slow_label_bound"],162)
        self.controls("aggregate_seed_chain","no_label_or_deck_factor_in_aggregate")

    def test_edge_series_bound_includes_boundary_and_interior(self):
        self.controls("edge_unweighting","edge_series_coefficient_controls_x_power")
        self.controls("edge_unweighting","interior_delta_one_is_covered")
        self.controls("edge_unweighting","exponential_series_term")

    def test_explicit_edge_constants_keep_width_and_fractional_order(self):
        self.assertEqual(edge_constant(0,Q(1,8)),1)
        self.assertEqual(edge_constant(2,Q(1,8)),128)
        self.assertEqual(edge_constant(3,Q(1,8)),32768)
        self.assertEqual(edge_constant(Q(3,2),Q(1,8)),128)
        self.controls("edge_unweighting","width_dependence_cannot_be_declared_uniform")

    def test_slow_and_auxiliary_normalization_factors_are_compatible(self):
        self.controls("matched_normalization","coefficient_conversion_exponents")
        self.controls("matched_normalization","radial_U_multiplier_bound")
        self.controls("matched_normalization","raw_slow_chart_multipliers")
        self.controls("matched_normalization","actual_band_window_slow_scale_ratio")

    def test_exact_amplitude_and_radial_chart_powers(self):
        self.controls("matched_normalization","amplitude_normalization_bound_eight")
        self.controls("matched_normalization","radial_U_normalization_bound_thirty_two")
        self.controls("matched_normalization","coarsest_index_is_common_minimum")

    def test_literal_native_carrier_bound_can_fail_for_finer_reference(self):
        self.controls("matched_normalization","finer_reference_literal_carrier_bound_fails")
        self.controls("matched_normalization","native_rounding_counterexample_is_exact")

    def test_time_floor_cost_changes_when_reference_and_index_are_mixed(self):
        self.controls("matched_normalization","native_time_floor_cannot_be_reused")
        self.controls("matched_normalization","arbitrary_time_floor_constant_suffices")
        self.controls("matched_normalization","rational_time_floor_majorants")

    def test_rescaling_preserves_positive_polynomial_coefficients(self):
        self.assertEqual(rescale(P({0:2,2:3}),Q(25,16)).record(),{"0":"2","2":"1875/256"})

    def test_coarsest_aggregate_matches_parent_covariance_envelope(self):
        r=aggregate(0,P(1),P(1),0,0,1)
        self.assertEqual(r["a"].evaluate(1),8)
        self.assertEqual(r["b"].evaluate(1),8)
        self.assertEqual(r["Omega"].evaluate(1),64)
        self.assertEqual((r["carrier_ratio"],r["time_floor"]),(1,6))

    def test_arbitrary_reference_keeps_additional_carrier_and_time_constants(self):
        r=aggregate(0,P(1),P(1),0,0,1,True)
        self.assertEqual(r["Omega"].evaluate(1),96)
        self.assertEqual((r["carrier_ratio"],r["time_floor"]),(2,1296))
        self.controls("aggregate_seed_chain","arbitrary_carrier_ratio_retained")

    def test_aggregate_bounds_substitute_into_all_initial_seed_scalars(self):
        self.controls("aggregate_seed_chain","aggregate_to_seed_substitution")

    def test_whole_pressure_moment_substitution_and_jet_counts(self):
        self.controls("aggregate_seed_chain","whole_pressure_moment_substitution")
        self.controls("aggregate_seed_chain","scalar_jet_ledger_unchanged")
        self.controls("aggregate_seed_chain","axial_moment_jet_ledger_unchanged")

    def test_invalid_edge_and_reference_data_rejected(self):
        for b,width in ((-1,1),(1,0),(1,2)):
            with self.assertRaises(ValueError):edge_constant(b,width)
        with self.assertRaises(TypeError):edge_constant(.5,1)
        with self.assertRaises(ValueError):rescale(P(1),0)
        with self.assertRaises(TypeError):rescale(1,1)
        with self.assertRaises(ValueError):aggregate(-1,P(1),P(1),0,0,1)
        with self.assertRaises(TypeError):aggregate(0,P(1),P(1),0,0,1,1)

    def test_pressure_memory_and_complete_covariance_are_retained(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))
        self.controls("complete_covariance","all_six_covariances_nonzero")
        self.controls("complete_covariance","omitted_axial_velocity_detected")

    def test_historical_bindings_and_frozen_record_reproduce(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_aggregate_join_not_promoted_to_per_label_or_residual_proof(self):
        self.assertTrue(self.record["flags"]["actual_atlas_geometry_not_presented_as_new"])
        for key in ("per_label_canonical_ODE_root_and_edge_estimates_independently_reproved","full_signed_sourced_seed_or_pointwise_residual_verified","uniform_infinite_stage_or_all_order_constant_proved","actual_profile_radius_tilt_width_or_accepted_tail_changed","unconditional_pulse_or_blowup_proved"):
            self.assertFalse(self.record["flags"][key])

    def test_scope_promotion_is_rejected(self):
        changed=deepcopy(self.record);changed["flags"]["unconditional_pulse_or_blowup_proved"]=True
        with self.assertRaises(ValueError):audit.validate(changed)

    def test_existing_frozen_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"evidence.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
