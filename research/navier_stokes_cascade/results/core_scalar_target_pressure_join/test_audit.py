#!/usr/bin/env python3
"""Current scalar jets, projected-source controls and complete pressure checks."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import review


class ScalarTargetJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_current_complete_angular_equation_moment(self):
        self.controls("current_scalar","current_scalar_is_integrated_complete_angular_equation")

    def test_Haar_divergence_and_both_velocity_moments(self):
        self.controls("current_scalar","current_Haar_mean_exact_divergence")
        self.controls("current_scalar","current_base_exact_divergence")
        self.controls("current_scalar","both_current_velocity_moments_zero")

    def test_projection_and_scalar_reconstruct_the_current_equation(self):
        for key in ("projector_weighted_mass_one","projected_source_zero_weighted_mass","projected_target_compact","two_channels_reconstruct_current_equation","both_current_channels_nonzero"):
            self.controls("current_scalar",key)

    def test_current_Haar_products_and_full_covariance_are_retained(self):
        for key in ("averaging_products_after_multiplication_required","complete_covariance_scalar_flux_retained","complete_covariance_radial_flux_retained"):
            self.controls("current_scalar",key)

    def test_pure_T_scalar_bound_uses_current_field_data(self):
        self.controls("current_scalar","pure_T_current_flux_bound")

    def test_mixed_scalar_norm_split_and_current_data_bound(self):
        self.controls("current_scalar","mixed_scalar_jet_split")
        self.controls("current_scalar","current_data_scalar_jet_bound")

    def test_all_mixed_scalar_jet_identities(self):
        self.controls("current_scalar","exact_scalar_jet")

    def test_third_and_fourth_product_jet_routing(self):
        self.controls("pressure_routing","exact_product_scalar_moment_routing")
        self.controls("pressure_routing","whole_K_dependence_retained")

    def test_both_pressure_mismatch_norm_budgets(self):
        self.controls("pressure_routing","axial_viscous_current_data_budget")
        self.controls("pressure_routing","slow_time_current_data_budget")

    def test_zero_scalar_derivative_does_not_erase_variable_K(self):
        self.controls("pressure_routing","scalar_M_zero_does_not_erase_K_pressure_response")

    def test_projection_alone_cannot_bound_scalar_moment(self):
        self.controls("projection_control","projected_target_zero_with_scalar")
        self.controls("projection_control","moving_projection_bump_retained")

    def test_zero_Haar_graph_curls_keep_both_moments_and_divergence(self):
        for key in ("old_and_fast_Haar_means_all_zero","old_and_fast_graph_divergence_exact","old_fast_new_both_Haar_moments_zero","graph_sine_curl_component_retained"):
            self.controls("fast_scalar",key)

    def test_postfast_scalar_has_old_new_crosses_and_fourth_jet(self):
        for key in ("fast_scalar_update_exact_with_all_three_products","fast_correlation_energy_positive","zero_Haar_does_not_preserve_scalar_defect","fourth_Z_scalar_jet_survives","old_new_cross_cannot_be_dropped"):
            self.controls("fast_scalar",key)

    def test_scalar_and_integrated_residual_have_different_epsilon_orders(self):
        self.controls("fast_scalar","no_extra_epsilon_for_scalar_change")
        self.controls("fast_scalar","moment_response_keeps_one_epsilon")

    def test_both_full_pressure_tails_and_complete_covariance(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))
        self.controls("complete_covariance","all_six_covariances_nonzero")
        self.controls("complete_covariance","omitted_axial_velocity_detected")

    def test_historical_bindings_match(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)

    def test_scientific_limits_are_preserved(self):
        for key in ("projected_target_bound_alone_controls_discarded_scalar_channel","operator_projection_control_is_actual_canonical_counterexample","actual_incoming_norms_or_complete_cycle_inputs_verified","actual_profile_or_accepted_tail_changed","unconditional_pulse_residual_improvement_or_blowup_proved"):
            self.assertFalse(self.record["flags"][key])

    def test_scope_tampering_rejected(self):
        r=deepcopy(self.record);r["flags"]["actual_incoming_norms_or_complete_cycle_inputs_verified"]=True
        with self.assertRaises(ValueError):audit.validate(r)

    def test_frozen_record_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_output_overwrite_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"evidence.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
