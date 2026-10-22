#!/usr/bin/env python3
"""Current target primitive, whole scalar moment and weighted-input tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import P,TargetBudget,fixture,primitive_recipe,rational_upper
from controls import review


class SignedTargetEdgeBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_exact_flat_tail_primitive_differential_identities(self):
        self.controls("two_edge_tail","exact_tail_differential_identity")
        self.controls("two_edge_tail","kernel_and_exact_tail_agree")

    def test_three_shell_regions_and_fractional_poles(self):
        self.controls("two_edge_tail","collars_inside_shell")
        self.controls("two_edge_tail","kernel_degree_dominates")
        self.assertEqual(primitive_recipe(Q(3,2),4,1)["kernel_degree"],0)
        self.assertEqual(primitive_recipe(4,4,1)["kernel_degree"],1)

    def test_rational_upper_retains_width_and_exponential_costs(self):
        self.controls("two_edge_tail","rational_exponential_envelope")
        self.controls("two_edge_tail","narrower_width_cost_retained")

    def test_giant_selected_constants_are_not_materialized_as_diagnostics(self):
        recipe=primitive_recipe(9,4,Q(1,10**30))
        self.assertEqual(recipe["width"],str(Q(1,10**30)))
        with self.assertRaises(ValueError):rational_upper(primitive_recipe(1,Q(1,1000),1),max_power=20)

    def test_current_full_equations_keep_divergence_and_both_velocity_moments(self):
        for fragment in ("current_base_divergence","current_mean_divergence","both_velocity_moments_zero"):
            self.controls("current_projection",fragment)

    def test_current_two_distinct_scalar_channels_are_exact(self):
        self.controls("current_projection","current_angular_moment_whole_scalar")
        self.controls("current_projection","current_axial_moment_whole_pressure_scalar")

    def test_whole_moving_pressure_product_cannot_be_frozen(self):
        for fragment in ("whole_pressure_compact_both_ends","moving_pressure_product_term_nonzero","omitted_moving_pressure_factor_detected"):
            self.controls("current_projection",fragment)

    def test_current_covariance_and_nonzero_old_haar_means_retained(self):
        for fragment in ("all_six_complete_covariances_retained","old_haar_velocities_not_reset","products_averaged_after_multiplication"):
            self.controls("current_projection",fragment)

    def test_both_moving_projectors_have_exact_weighted_normalization(self):
        self.controls("current_projection","moving_projector_has_unit_weighted_mass")
        self.controls("current_projection","projected_source_exact_zero_moment")

    def test_haar_shifted_map_is_ordinary_only_after_exact_projection(self):
        for fragment in ("Haar_shifted_primitive_equals_ordinary","signed_target_compact_both_ends","signed_target_reconstructs_source"):
            self.controls("current_projection",fragment)

    def test_missing_moment_projection_loses_outer_compactness(self):
        self.controls("current_projection","omitted_projection_loses_outer_compactness")

    def test_scaled_two_bump_source_has_exact_zero_total_and_nonzero_primitive(self):
        for fragment in ("scaled_bump_mass_is_positive","observation_between_two_disjoint_bumps","exact_two_bump_total_zero","first_bump_primitive_mass"):
            self.controls("unweighted_obstruction",fragment)

    def test_exact_exponential_gap_defeats_any_fixed_edge_pole(self):
        self.controls("unweighted_obstruction","exponential_weight_gap")
        self.controls("unweighted_obstruction","any_fixed_weighted_pole_is_defeated")

    def test_small_all_fixed_unweighted_jets_do_not_supply_missing_weight(self):
        self.controls("unweighted_obstruction","all_fixed_unweighted_derivatives_can_tend_zero")
        self.controls("unweighted_obstruction","exact_moment_not_promoted_to_weighted_bound")

    def test_projection_and_two_edge_positive_polynomial_substitution(self):
        self.controls("target_root_chain","current_projection_budget")
        self.controls("target_root_chain","two_edge_primitive_budget")

    def test_reciprocal_weight_pole_fits_fixed_order_root_input(self):
        self.controls("target_root_chain","normalized_target_fits_root_input")
        self.assertEqual(fixture().edge(5)["d"],18)

    def test_nonlinear_chart_composition_keeps_higher_mixed_jets(self):
        self.controls("chart_composition","nonlinear_mixed_chart_composition")
        self.controls("chart_composition","dropping_higher_chart_derivatives_detected")

    def test_no_extra_source_jet_in_primitive_or_chart_conversion(self):
        self.controls("target_root_chain","primitive_and_weight_division_add_no_source_jet")

    def test_signed_curl_scalar_channel_needs_underlying_n_plus_13(self):
        self.controls("target_root_chain","signed_curl_moment_channel_uses_n_plus_13")
        self.controls("target_root_chain","scalar_moment_channel_is_not_full_row_proof")
        self.controls("target_root_chain","initial_moment_epsilon_three_is_stronger_than_signed_alpha")

    def test_original_pressure_memory_stays_in_replayed_controls(self):
        self.assertTrue(all(self.result["parts"]["pressure_memory"]["checks"].values()))

    def test_invalid_shell_width_or_pole_inputs_rejected(self):
        for p,L,w in ((-1,4,1),(1,0,1),(1,4,0),(1,4,2)):
            with self.assertRaises(ValueError):primitive_recipe(p,L,w)
        with self.assertRaises(TypeError):primitive_recipe(.5,4,1)
        row=fixture();row.poles=lambda n:-1
        with self.assertRaises(ValueError):row.edge(2)
        row.poles=lambda n:3-n
        with self.assertRaises(ValueError):row.edge(2)

    def test_nonpolynomial_bound_data_are_rejected(self):
        row=fixture();row.row=lambda n:1
        with self.assertRaises(TypeError):row.H(1)

    def test_historical_bindings_and_frozen_record_reproduce(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),5)
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_conditional_weighted_gate_not_promoted_to_actual_residual_or_blowup(self):
        for key,value in self.record["flags"].items():
            if not value:
                changed=deepcopy(self.record);changed["flags"][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)

    def test_missing_history_or_remaining_inequality_is_rejected(self):
        for field in ("source_hashes","remaining_obligations"):
            changed=deepcopy(self.record)
            if isinstance(changed[field],dict):changed[field].pop(next(iter(changed[field])))
            else:changed[field].pop(0)
            with self.assertRaises(ValueError):audit.validate(changed)

    def test_existing_evidence_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"keep.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
