"""Exact column integrals, all moving-map terms and unchanged canonical root tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,toy,original_fixture
from bounds import P,HaarColumnBudget,from_original,root_join,entrance

class CanonicalHaarColumnTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        vals=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(vals,key);self.assertTrue(all(vals),key)
    def test_all_six_unrooted_principal_tensor_entries_kept(self):self.control('all_six_principal_tensor_entries')
    def test_selected_columns_use_original_full_integrals(self):self.control('two_selected_columns');self.control('two_sign_columns_give_nonzero')
    def test_whole_interval_tensor_jets_through_order_three(self):self.control('all_tensor_entry_jets_enclosed')
    def test_matrix_infinity_norm_keeps_both_sign_columns(self):self.control('induced_matrix_row_sum_jet_bound')
    def test_actual_higher_column_jets_are_not_zero(self):self.control('higher_column_jets_are_nonzero')
    def test_temporal_cutoff_remains_in_exact_column(self):self.control('discarding_temporal_cutoff_changes')
    def test_Gaussian_derivatives_remain_in_pulse_jet_input(self):self.control('differentiating_only_t_over_P')
    def test_physical_radius_is_removed_only_after_s0_normalization(self):self.control('angular_half_transverse_mass');self.control('radius_independence_only_after')
    def test_angular_half_is_not_dropped(self):self.control('omitting_angular_half_doubles')
    def test_complete_moving_translation_Haar_derivatives(self):self.control('moving_translation_derivatives_cancel')
    def test_moving_terms_do_not_vanish_pointwise(self):self.control('moving_center_terms_are_nonzero')
    def test_translation_requires_complete_periodic_domain(self):self.control('truncating_to_nonperiodic_domain')
    def test_original_Jacobian_and_normalized_cover(self):self.control('original_Jacobian');self.control('normalized_Haar_needs_no_deck_count')
    def test_factorial_normalization_in_fixed_coefficient_chart(self):self.control('raw_multivariate_factorial_count')
    def test_original_factory_derives_curvature_from_actual_tilt(self):self.control('actual_fixed_tilt_not_replaced');self.control('original_Gaussian_curvature')
    def test_original_band_and_tangent_gates_are_not_reselected(self):self.control('original_normal_tangent_and_band_gates')
    def test_joint_H_and_W_inputs_use_same_canonical_pulse(self):self.control('joint_H_input');self.control('joint_W_input')
    def test_value_gaps_targets_and_full_curl_inputs_unchanged(self):self.control('all_canonical_value_gaps');self.control('target_curl_loghalf_and_trace')
    def test_join_supplies_inverse_root_signed_and_complete_curl_recipes(self):self.control('joined_inverse_root_signed_and_complete_curl')
    def test_Gaussian_square_integral_rational_majorant(self):self.control('rational_Gaussian_integral_upper')
    def test_no_uniform_curvature_from_tilt_upper_bound_alone(self):self.control('fixed_positive_tilt_required')
    def test_zero_or_invalid_actual_tilt_is_rejected(self):
        p,sp,ell,kw=original_fixture()
        for u in (0,-1,1):
            with self.assertRaises(ValueError):from_original(p,sp,ell,**dict(kw,tilt=u))
        with self.assertRaises(TypeError):from_original(p,sp,ell,**dict(kw,tilt=0.01))
    def test_unknown_early_band_is_not_reselected(self):
        p,sp,ell,kw=original_fixture()
        with self.assertRaises(ValueError):from_original(p,sp,16,**kw)
    def test_invalid_column_inputs_or_noncanonical_root_join_are_rejected(self):
        h=toy()[-1].homogeneous
        for c in (0,-1):
            with self.assertRaises(ValueError):HaarColumnBudget(h,c)
        with self.assertRaises(TypeError):HaarColumnBudget(h,0.5)
        for d in (0,2,3.5):
            with self.assertRaises(ValueError):HaarColumnBudget(h,1,d)
        with self.assertRaises(TypeError):HaarColumnBudget(object(),1)
        with self.assertRaises(TypeError):root_join(entrance.inverse.roots.fixture(),toy()[-1],lambda j:P(1))
    def test_invalid_jet_orders_are_rejected(self):
        c=toy()[-1]
        for j in (-1,Q(1,2)):
            with self.assertRaises((ValueError,TypeError)):c.matrix(j)
    def test_all_exact_controls_and_historical_source_hashes_pass(self):self.assertTrue(all(self.r['checks'].values()));self.assertTrue(all(self.record['gates'].values()))
    def test_eight_frozen_parent_records_preserved(self):self.assertEqual(len(self.record['dependency_evidence_sha256']),8);self.assertEqual(self.record['source_hashes'],audit.sources())
    def test_principal_tensor_not_promoted_to_complete_covariance_or_blowup(self):
        for name in ('canonical_and_inherited_families_or_principal_and_full_curl_covariances_identified','actual_full_background_suprema_weighted_target_or_full_curl_operator_input_certified','unconditional_pulse_complete_residual_infinite_field_or_blowup_proved'):
            self.assertIs(self.record['flags'][name],False)
        self.assertEqual(len(self.record['remaining_obligations']),4)
    def test_tampered_scope_recipe_and_source_provenance_rejected(self):
        for field in ('flags','review','source_hashes'):
            r=deepcopy(self.record);r[field]['tampered']=True
            with self.assertRaises(ValueError):audit.validate(r)
    def test_frozen_evidence_reproduces(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()))
    def test_cli_uses_fresh_output_and_rejects_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/'replay.json';audit.main(['--verify-record',str(Path(__file__).parent/'evidence.json'),'--output',str(dest)])
            self.assertEqual(json.loads(dest.read_text()),self.record)
            with self.assertRaises(SystemExit):audit.main(['--output',str(dest)])
if __name__=='__main__':unittest.main()
