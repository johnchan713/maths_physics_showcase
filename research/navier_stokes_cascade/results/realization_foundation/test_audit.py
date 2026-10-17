#!/usr/bin/env python3
"""Current representatives, whole-interval diagonal and exact cutoff controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from realized_bounds import B0,GAIN,RHO,KAPPA,orders,gain,tuple_loss,cutoff_loss,field_loss,comparison_budget,representative_review,scale_review
from realized_checks import implicit_review,polynomial_review,apply_word,canonical,add,mul,grad,ONE,T,Z,power
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main

class RealizationFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_parent_is_exact_live_foundation(self):self.assertEqual(protocol()['parent_commit'],'07b30739332662b1bca55018fa735bb0391ae5d5')
    def test_current_indexing_identity_for_entire_affine_sequence(self):self.assertEqual(B0-GAIN,RHO)
    def test_every_representative_has_same_positive_stage_slope(self):self.assertEqual(orders(2),{k:v+GAIN for k,v in orders(1).items()})
    def test_weakest_all_stage_margin_is_one_minus_two_kappa(self):self.assertTrue(representative_review()['checks']['weakest_margin_is_one_minus_two_kappa'])
    def test_direct_swirl_is_kept_separate_from_compact_potential(self):self.assertIn('direct_axisymmetric_swirl',orders(1));self.assertFalse(FLAGS['compact_full_direct_swirl_vector_potential_assumed'])
    def test_complete_pressure_after_all_four_operations_and_fixed_gauge(self):self.assertIn('fixed compact gauge',protocol()['representatives']);self.assertIn('after all four',protocol()['pressure'])
    def test_cutoff_loss_reserves_one_full_curl(self):self.assertEqual(cutoff_loss(3),tuple_loss(4))
    def test_finite_state_loss_reserves_fixed_log_absorption(self):self.assertEqual(field_loss(3),cutoff_loss(3)+1)
    def test_bad_stage_or_jet_indices_rejected(self):
        for n in (-1,True,F(1,2)):
            with self.assertRaises(ValueError):gain(n)
            with self.assertRaises(ValueError):tuple_loss(n)
        with self.assertRaises(ValueError):gain(0)
    def test_original_h_domain_retained(self):
        for h in (0,F(1,100),-1):
            with self.assertRaises(ValueError):comparison_budget(h=h)
    def test_mixed_actual_q_cutoff_operators_commute_in_all_226_coefficients(self):self.assertEqual(implicit_review()['mixed_coefficient_comparisons'],226)
    def test_all_54_positive_jets_use_positive_cutoff_derivatives(self):self.assertEqual(implicit_review()['positive_derivative_cases'],54)
    def test_exact_z_and_time_cutoff_signs(self):
        self.assertEqual(apply_word('z',F(1,200))[0],{1:{(1,1):F(2)}})
        self.assertEqual(apply_word('t',F(1,200))[0],{1:{(0,1):F(-1)}})
    def test_canonical_eta_denominator_identity(self):
        h=F(1,200);self.assertFalse(canonical({(0,0):F(1),(0,1):F(-1),(2,1):2*h},h))
    def test_full_cut_curl_exactly_divergence_free(self):self.assertTrue(polynomial_review()['checks']['complete_cut_velocity_exactly_divergence_free'])
    def test_naive_velocity_cutoff_failure_detected(self):self.assertTrue(polynomial_review()['checks']['naive_velocity_cutoff_breaks_divergence'])
    def test_axisymmetric_direct_swirl_divergence_free(self):self.assertTrue(polynomial_review()['checks']['axisymmetric_direct_swirl_exactly_divergence_free'])
    def test_pressure_and_viscous_cutoff_terms_nonzero(self):
        c=polynomial_review()['checks'];self.assertTrue(c['pressure_cutoff_gradient_is_nonzero']);self.assertTrue(c['viscous_cutoff_commutators_are_nonzero'])
    def test_time_pressure_gauge_cannot_be_changed_before_cutoff(self):
        gauge=add(ONE,T);cut=add(ONE,power(Z,2));self.assertFalse(any(grad(gauge)));self.assertTrue(any(grad(mul(cut,gauge))))
    def test_all_six_residual_groups_are_individually_necessary(self):self.assertEqual(polynomial_review()['nonzero_individually_necessary_terms'],['new_old','new_square','old_new','pressure','time','viscosity'])
    def test_polynomial_fixture_is_not_actual_supported_giant_field(self):self.assertFalse(polynomial_review()['actual_supported_giant_field_evaluated']);self.assertTrue(polynomial_review()['local_polynomial_cutoff_not_claimed_C_infinity_compact'])
    def test_whole_interval_scales_cover_extreme_synthetic_constants(self):self.assertTrue(all(scale_review()['checks'].values()))
    def test_fixed_stage_comparison_all_targets_and_jet_orders(self):
        for m,N,h in ((0,0,F(1,200)),(1,3,F(9,1000)),(20,2,F(1,10000))):self.assertTrue(all(comparison_budget(m,N,h)['checks'].values()))
    def test_coarse_zero_jet_budget_is_12816_unevaluated_stages(self):
        r=comparison_budget();self.assertEqual(r['comparison_stage'],12816);self.assertFalse(r['actual_fields_or_scales_evaluated'])
    def test_no_fixed_q_or_arbitrary_j_of_q_claim(self):self.assertFalse(FLAGS['arbitrary_j_of_q_substitution_justified']);self.assertFalse(FLAGS['uniform_stage_constants_or_raw_uncut_series_convergence_proved'])
    def test_all_961_historical_and_18_direct_source_bindings(self):self.assertEqual(len(self.record['inherited_provenance']),961);self.assertEqual(len(self.record['source_hashes']),18)
    def test_frozen_record_and_roundtrip_match(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r);validate(json.loads(json.dumps(self.record)))
    def test_no_endpoint_global_or_unconditional_scope_promotion(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_cached_exact_reviews_not_mutated_by_generated_record(self):
        r=generate();r['reviews']['actual_implicit_q_all_order_cutoff_mechanism']['checks']['all_mixed_cutoff_derivatives_commute_exactly']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(implicit_review()['checks']['all_mixed_cutoff_derivatives_commute_exactly'])
    def test_same_count_provenance_substitution_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])
    def test_protocol_gates_sources_remaining_cannot_be_erased(self):
        for k in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[k],dict):r[k].pop(next(iter(r[k])))
            else:r[k].pop()
            with self.assertRaises(ValueError):validate(r)
    def test_cli_never_overwrites_existing_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')

if __name__=='__main__':unittest.main()
