"""Full cylindrical curl, coefficient jets, chart transport and scope regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,manufactured,original_fixture
from bounds import P,CurlBudget,from_original,join_root,entrance

class CompleteCurlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        vals=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(vals,key);self.assertTrue(all(vals),key)
    def test_full_cylindrical_graph_operator_identity(self):self.control('full_cylindrical_graph_coefficient_identity');self.control('independent_component_curl')
    def test_complete_graph_curl_is_divergence_free(self):self.control('full_graph_curl_is_divergence_free')
    def test_transversality_uses_full_normal(self):self.control('transverse_input_uses_full_normal')
    def test_cylindrical_connection_not_deleted(self):self.control('cylindrical_connection_is_nonzero')
    def test_radial_graph_terms_not_deleted(self):self.control('radial_graph_terms_are_nonzero')
    def test_slow_axial_terms_not_deleted(self):self.control('slow_axial_terms_are_nonzero')
    def test_differentiated_normal_potential_retained(self):self.control('differentiated_normal_potential');self.control('normal_reciprocal_derivatives')
    def test_graph_coefficient_derivatives_retained(self):self.control('graph_coefficient_derivative_is_not_zero')
    def test_whole_box_sum_of_every_matrix_jet(self):self.control('sum_of_all_matrix_raw_jets_enclosed')
    def test_conjugated_coefficient_jets_on_same_box(self):self.control('conjugated_U_coefficient_raw_jets')
    def test_D0_retains_extra_normal_jet(self):self.control('D0_uses_one_additional')
    def test_scale_normalization_has_no_remaining_negative_epsilon(self):self.control('native_scale_normalization_keeps_graph')
    def test_exact_U_operator_conjugation(self):self.control('U_operator_requires_conjugation')
    def test_wrong_unchanged_R_matrix_in_U_is_detected(self):self.control('reusing_R_first_order_matrix_in_U')
    def test_finite_composition_matches_independent_Taylor(self):self.control('finite_composition_matches_direct_Taylor');self.control('positive_composition_majorant')
    def test_inverse_power_falling_factorials_bound_all_test_orders(self):self.control('falling_factor_bound')
    def test_original_exponent_bounds_are_sufficient(self):self.control('original_exponent_interval')
    def test_production_factory_keeps_original_band_and_gates(self):self.control('original_same_band_normal_tangent')
    def test_H_W_and_curl_share_full_normal_recipe(self):self.control('potential_lift_is_the_same')
    def test_same_R_root_recipes_and_consistent_U_transport(self):self.control('same_R_recipe_and_all_U_jets_transported');self.control('complete_curl_callback_derived')
    def test_root_value_gaps_mass_radius_and_width_unchanged(self):self.control('value_gaps_mass_width_radius_retained')
    def test_complete_curl_and_signed_maps_use_joined_coefficients(self):self.control('complete_curl_and_signed_root_recipes')
    def test_external_signed_target_is_explicitly_transported(self):self.control('separately_supplied_signed_target_is_also_transported')
    def test_full_axial_pressure_moment_background_order_explicit(self):self.control('whole_moving_axial_moment_retains_n_plus_thirteen')
    def test_invalid_hulls_inputs_or_mixed_chart_names_rejected(self):
        h=manufactured()[-1].homogeneous
        for args in ((0,1),(2,1),(-1,1)):
            with self.assertRaises(ValueError):CurlBudget(h,*args)
        with self.assertRaises(TypeError):CurlBudget(h,0.5,1)
        with self.assertRaises(TypeError):CurlBudget(object(),1,2)
        with self.assertRaises(ValueError):CurlBudget(h,1,2,dimension=4)
        op=CurlBudget(h,1,2)
        with self.assertRaises(ValueError):op.normalized(0,'mixed')
        with self.assertRaises(ValueError):op.transport(lambda j:P(1),'mixed')
        with self.assertRaises(ValueError):join_root(entrance.inverse.roots.fixture(),op,'mixed')
        with self.assertRaises(TypeError):join_root(object(),op)
    def test_invalid_derivative_orders_rejected(self):
        op=manufactured()[-1]
        for n in (-1,Q(1,2)):
            with self.assertRaises((ValueError,TypeError)):op.normalized(n)
        with self.assertRaises(ValueError):op.inverse_map(0)
    def test_unknown_early_band_not_reselected(self):
        p,sp,ell,kw=original_fixture()
        with self.assertRaises(ValueError):from_original(p,sp,16,**kw)
    def test_all_exact_controls_and_historical_sources_pass(self):self.assertTrue(all(self.r['checks'].values()));self.assertTrue(all(self.record['gates'].values()))
    def test_eight_frozen_dependencies_and_current_sources(self):self.assertEqual(len(self.record['dependency_evidence_sha256']),8);self.assertEqual(self.record['source_hashes'],audit.sources())
    def test_no_actual_background_or_residual_or_blowup_promotion(self):
        for k in ('actual_background_suprema_weighted_target_or_cutoff_input_certified','unconditional_pulse_complete_residual_infinite_field_or_blowup_proved'):
            self.assertIs(self.record['flags'][k],False)
        self.assertTrue(self.record['flags']['all_six_covariances_axial_mean_and_reserved_pressure_support_obligations_retained'])
    def test_tampered_scope_recipe_or_hashes_rejected(self):
        for field in ('flags','review','source_hashes'):
            r=deepcopy(self.record);r[field]['tampered']=True
            with self.assertRaises(ValueError):audit.validate(r)
    def test_frozen_evidence_reproduces(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()))
    def test_cli_replays_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/'replay.json';audit.main(['--verify-record',str(Path(__file__).parent/'evidence.json'),'--output',str(dest)])
            self.assertEqual(json.loads(dest.read_text()),self.record)
            with self.assertRaises(SystemExit):audit.main(['--output',str(dest)])
if __name__=='__main__':unittest.main()
