#!/usr/bin/env python3
"""Complete initial Haar rows, weighted maps and rate/scope regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import P,KEYS,row_budget,weighted_covariance,fast_weighted,weighted_operator_budget,fixture,root
from controls import review


class InitialHaarRowBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_zero_Haar_fast_and_nonzero_slow_state_kept_distinct(self):
        self.controls("current_complete_equations","actual_fast_all_Haar_means_zero")
        self.controls("current_complete_equations","actual_slow_Haar_means_nonzero")

    def test_coefficient_divergence_and_both_scalar_velocity_moments(self):
        self.controls("current_complete_equations","fast_and_slow_divergence_exact")
        self.controls("current_complete_equations","both_velocity_moments_zero")

    def test_full_radial_source_and_current_pressure_reconstructed(self):
        for fragment in ("full_radial_equation_matches_split","full_pressure_defect_matches_split","full_pressure_profile_matches_split"):
            self.controls("current_complete_equations",fragment)

    def test_both_full_current_Haar_tangential_rows_match_independent_split(self):
        self.controls("current_complete_equations","full_current_angular_row_matches_split")
        self.controls("current_complete_equations","full_current_axial_row_matches_split")

    def test_fast_linear_terms_vanish_only_in_Haar_rows(self):
        self.controls("current_complete_equations","fast_linear_time_and_viscosity_Haar_zero")
        self.controls("sharper_Haar_rate","unaveraged_fast_graph_viscosity_not_promoted")

    def test_fast_slow_products_are_not_deleted_pointwise(self):
        self.controls("current_complete_equations","fast_slow_crosses_zero_only_after_averaging")

    def test_complete_covariance_and_base_radial_velocity_remain(self):
        self.controls("current_complete_equations","all_six_complete_covariances_nonzero")
        self.controls("current_complete_equations","base_radial_velocity_not_zero")

    def test_all_fast_quadratic_correlations_matter(self):
        for fragment in ("omitted_fast_pressure_squares_detected","omitted_fast_angular_correlations_detected","omitted_fast_axial_correlations_detected"):
            self.controls("current_complete_equations",fragment)

    def test_current_compact_pressure_and_two_whole_scalar_identities(self):
        for fragment in ("current_pressure_compact_both_ends","angular_scalar_identity_retained","whole_axial_pressure_scalar_identity_retained","moving_pressure_factor_cannot_be_frozen"):
            self.controls("current_complete_equations",fragment)

    def test_nonzero_old_Haar_state_cannot_use_initial_zero_Haar_split(self):
        self.controls("current_complete_equations","nonzero_old_Haar_breaks_zero_Haar_split")

    def test_all_retained_groups_fit_sharper_Haar_exponent(self):
        self.controls("sharper_Haar_rate","fits_sharper_Haar_rate")

    def test_pressure_itself_not_promoted_to_improved_row_order(self):
        self.controls("sharper_Haar_rate","pressure_itself_has_order_one")
        self.controls("sharper_Haar_rate","full_mean_pressure_fast_square_is_higher")

    def test_explicit_polynomials_reproduce_exact_numeric_substitution(self):
        self.controls("weighted_norm_chain","full_positive_row_substitution")

    def test_norm_budgets_retain_fast_correlations_and_whole_pressure(self):
        self.controls("weighted_norm_chain","fast_correlations_in_both_row_budgets")
        self.controls("weighted_norm_chain","whole_pressure_not_removed_from_axial_budget")

    def test_current_pressure_uses_n_plus_one_source_not_unavailable_N_jets(self):
        self.controls("weighted_norm_chain","pressure_jet_input_is_n_plus_one_not_N_plus_two")

    def test_five_variable_count_and_actual_fast_A1_are_retained(self):
        self.controls("weighted_norm_chain","five_variable_raw_to_normalized_count")
        self.controls("weighted_norm_chain","actual_axial_A1_budget_retained")
        self.controls("weighted_norm_chain","weighted_fast_radial_extra_input_jet")

    def test_weighted_covariance_not_replaced_by_edge_unweighted_bound(self):
        self.controls("weighted_norm_chain","weighted_primary_norms_not_unweighted")
        self.controls("weighted_norm_chain","real_selected_curl_square_is_not_full_covariance")

    def test_derived_current_rows_feed_existing_two_edge_original_root_chain(self):
        self.controls("weighted_norm_chain","current_rows_feed_two_edge_and_original_root_chain")

    def test_edge_directed_compact_map_handles_nonzero_moments_without_erasing_Ae(self):
        for fragment in ("nonzero_moment_left_primitive_has_wrong_outer_end","edge_directed_cutoff_primitive_is_compact","interior_compact_error_retained"):
            self.controls("weighted_fixed_maps",fragment)

    def test_weighted_map_constants_retain_interior_gap_and_weight_floor(self):
        self.controls("weighted_fixed_maps","weighted_compact_constants_retain_profile_and_gap")
        self.controls("weighted_fixed_maps","weighted_A1_not_assigned_unweighted_constant")
        self.controls("weighted_fixed_maps","primitive_collar_must_fit_actual_cutoff_plateau")
        with self.assertRaises(ValueError):weighted_operator_budget(1,4,1,P(1),P(1),0)

    def test_both_exceptional_pressure_tails_remain(self):
        self.assertTrue(all(self.result["parts"]["complete_pressure_tails"]["checks"].values()))

    def test_missing_negative_or_mixed_bound_inputs_rejected(self):
        with self.assertRaises(ValueError):row_budget(0,{})
        with self.assertRaises(ValueError):row_budget(-1,fixture())
        data={k:Q(1) for k in KEYS};data["b"]=-1
        with self.assertRaises(ValueError):row_budget(0,data)
        data=fixture();data["b"]=1
        with self.assertRaises(TypeError):row_budget(0,data)

    def test_six_frozen_dependencies_and_record_reproduce(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_no_unaveraged_later_state_residual_or_blowup_promotion(self):
        for key,value in self.record["flags"].items():
            if not value:
                changed=deepcopy(self.record);changed["flags"][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)

    def test_missing_history_or_remaining_obligations_rejected(self):
        for field in ("source_hashes","remaining_obligations"):
            changed=deepcopy(self.record)
            if isinstance(changed[field],dict):changed[field].pop(next(iter(changed[field])))
            else:changed[field].pop(0)
            with self.assertRaises(ValueError):audit.validate(changed)

    def test_existing_evidence_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"keep.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
