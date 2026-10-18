#!/usr/bin/env python3
"""Exact canonical matching, diagonal majorant and missing-gauge controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
import audit
from realization import B0,DELTA,RHO,K,index,parameters,order,gain,orders,tuple_loss,cut_loss,total_loss,comparison,choose_log_scale,scale_review,universal_controls
from gauge_checks import review as gauge_review

class CanonicalRealizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate()
    def test_current_canonical_parent(self):self.assertEqual(audit.protocol()['parent_commit'],'20bb5a62a4770cc93fb7372705b1991d03905a4e')
    def test_canonical_seed_and_gain_match_own_parent(self):self.assertTrue(all(self.record['reviews']['canonical_finite_input']['checks'].values()))
    def test_stage_identity_for_all_affine_stages(self):self.assertEqual(B0-DELTA,RHO)
    def test_every_representative_exact_slope(self):self.assertEqual(orders(2),{k:v+DELTA for k,v in orders(1).items()})
    def test_weakest_representative_margin(self):self.assertEqual(min(orders(1).values())-DELTA,1-2*K)
    def test_direct_swirl_has_separate_representative(self):self.assertIn('direct_axisymmetric_swirl',orders(1))
    def test_complete_pressure_difference_in_ledger(self):self.assertIn('complete_fixed_gauge_pressure_difference',orders(1))
    def test_extra_curl_loss_reserved(self):self.assertEqual(cut_loss(3),tuple_loss(4))
    def test_total_field_loss_absorbs_fixed_logs(self):self.assertEqual(total_loss(3),cut_loss(3)+1)
    def test_invalid_indices_rejected(self):
        for x in (-1,True,F(1,2),1.0):
            with self.assertRaises(ValueError):index(x)
    def test_original_h_domain_retained(self):
        for h in (0,-1,F(1,100),F(1,50)):
            with self.assertRaises(ValueError):parameters(h)
    def test_whole_interval_rational_scale_majorants(self):self.assertTrue(all(scale_review()['checks'].values()))
    def test_scale_exact_fraction_for_extreme_constants(self):self.assertIsInstance(choose_log_scale(1,F(1,10**9),[(10**60,10**30)]),F)
    def test_scale_threshold_is_strict_not_rounded_endpoint(self):
        g=F(1,300);s=choose_log_scale(7,g,[(13,40)],100,20)
        self.assertLess(13+F(40)**2/(2*g)-g*s/2+7,0);self.assertGreater(s,101)
    def test_empty_scale_constraints_rejected(self):
        with self.assertRaises(ValueError):choose_log_scale(1,1,[])
    def test_nonpositive_gamma_rejected(self):
        for gamma in (0,-1):
            with self.assertRaises(ValueError):choose_log_scale(1,gamma,[(1,2)])
    def test_negative_scale_inputs_rejected(self):
        for pairs,previous,domain in (([(-1,0)],0,2),([(1,-2)],0,2),([(0,0)],-1,2),([(0,0)],0,-1)):
            with self.assertRaises(ValueError):choose_log_scale(1,1,pairs,previous,domain)
    def test_fixed_stage_budget_zero_jet(self):self.assertEqual(comparison()['comparison_stage'],12816)
    def test_fixed_stage_budget_all_six_residual_groups(self):
        for m,N,h in ((0,0,F(1,200)),(5,11,F(9,1000)),(20,2,F(1,10000))):self.assertTrue(all(comparison(m,N,h)['checks'].values()))
    def test_comparison_not_actual_field_evaluation(self):self.assertFalse(self.record['reviews']['fixed_stage_complete_residual_comparison']['actual_fields_or_scales_evaluated'])
    def test_actual_implicit_q_mixed_coefficients(self):self.assertEqual(universal_controls().implicit_review()['mixed_coefficient_comparisons'],226)
    def test_actual_q_positive_cutoff_jets(self):self.assertEqual(universal_controls().implicit_review()['positive_derivative_cases'],54)
    def test_complete_cut_curl_divergence(self):self.assertTrue(universal_controls().polynomial_review()['checks']['complete_cut_velocity_exactly_divergence_free'])
    def test_all_six_necessary_residual_groups(self):self.assertEqual(universal_controls().polynomial_review()['nonzero_individually_necessary_terms'],['new_old','new_square','old_new','pressure','time','viscosity'])
    def test_uncut_potential_gauge_invisible(self):self.assertTrue(gauge_review()['checks']['uncut_gradient_gauge_has_zero_curl'])
    def test_cut_potential_gauge_artifact_nonzero(self):self.assertGreater(gauge_review()['potential_gauge_artifact_components'],0)
    def test_potential_gauge_artifact_exact_product_rule(self):self.assertTrue(gauge_review()['checks']['exact_potential_gauge_artifact_is_grad_chi_cross_grad_phi'])
    def test_gauge_artifact_can_be_divergence_free_and_nonzero(self):self.assertTrue(gauge_review()['checks']['potential_gauge_artifact_still_divergence_free']);self.assertTrue(gauge_review()['checks']['potential_gauge_after_cutoff_changes_velocity'])
    def test_common_time_pressure_gauge_harmless(self):self.assertTrue(gauge_review()['checks']['common_time_gauge_cancels_from_all_pressure_increments'])
    def test_actual_q_gauge_counterexample_has_negative_power(self):
        r=gauge_review()['actual_q_gauge_counterexample'];self.assertEqual(F(r['pressure_gradient_factor_per_chi_prime']),F(100,133));self.assertEqual(F(r['power']),-F(99,200));self.assertFalse(r['actual_canonical_fields'])
    def test_independent_time_pressure_gauges_not_harmless(self):self.assertGreater(gauge_review()['stage_pressure_gauge_artifact_components'],0)
    def test_premature_middle_pressure_truncation_detected(self):self.assertTrue(gauge_review()['checks']['prematurely_discarding_middle_pressure_change_is_detected'])
    def test_controls_not_actual_compact_fields(self):self.assertFalse(gauge_review()['actual_supported_gauge_sequence_evaluated'])
    def test_frozen_record_exact_roundtrip(self):
        f=json.loads((audit.HERE/'evidence.json').read_text());audit.validate(f);audit.compare_records(self.record,f);audit.compare_records(self.record,json.loads(json.dumps(self.record)))
    def test_old_realization_not_identified(self):self.assertFalse(audit.FLAGS['old_realized_field_identified_with_canonical_field'])
    def test_no_endpoint_or_unconditional_promotion(self):
        for flag in ('canonical_closed_time_endpoint_atlas_joined','canonical_global_force_energy_and_comparison_joined','unconditional_complete_residual_improvement_or_blowup_proved','formal_certificate'):
            bad=deepcopy(self.record);bad['flags'][flag]=True
            with self.assertRaises(ValueError):audit.validate(bad)
    def test_source_hash_tampering_rejected(self):
        bad=deepcopy(self.record);key=next(iter(bad['source_hashes']));bad['source_hashes'][key]='0'*64
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_same_count_history_key_substitution_rejected(self):
        bad=deepcopy(self.record);v=bad['inherited_provenance']['bindings'];v['forged']=v.pop(next(iter(v)))
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_scale_majorant_tampering_rejected(self):
        bad=deepcopy(self.record);bad['reviews']['rational_whole_interval_scale']['rows'][0]['whole_interval_log_margin_upper']='1'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_generated_reviews_do_not_mutate_cache(self):
        bad=deepcopy(self.record);bad['reviews']['potential_and_pressure_gauge_controls']['checks'].clear()
        self.assertTrue(audit.reviews()['potential_and_pressure_gauge_controls']['checks'])
    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'kept.json';p.write_text('keep')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
if __name__=='__main__':unittest.main()
