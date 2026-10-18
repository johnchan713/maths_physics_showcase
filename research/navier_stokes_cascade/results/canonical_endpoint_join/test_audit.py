#!/usr/bin/env python3
"""Canonical trace matching, sharp factorial bounds and exact obstruction."""
from copy import deepcopy
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import tempfile
import unittest
import audit
from endpoint_bounds import universal_bounds,universal_checks,gevrey_row,gevrey_review
from endpoint_checks import dt,add,scale,obstruction_review

class CanonicalEndpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate()
    def test_same_canonical_local_parent(self):self.assertEqual(audit.protocol()['parent_commit'],'4b9bdd40944f851fc0aebc10fb47bf22f30b4e37')
    def test_canonical_own_entrance_roots_and_current_sources(self):self.assertTrue(audit.FLAGS['own_canonical_primary_entrance_roots_and_current_generated_sources_used'])
    def test_no_old_realized_field_identification(self):self.assertFalse(audit.FLAGS['old_endpoint_or_global_field_identified_with_canonical_field'])
    def test_all_canonical_input_checks(self):self.assertTrue(all(self.record['reviews']['canonical_current_realization_input']['checks'].values()))
    def test_closed_jet_join_not_velocity_only(self):self.assertIn('including NEXT',audit.protocol()['blocking_inequality'])
    def test_enlarged_full_boundary_window(self):self.assertEqual(universal_bounds().window(F(2)**-24,F(2)**-20),dict(first=18,last=26,count=9))
    def test_old_narrow_window_not_reused(self):self.assertEqual(universal_bounds().window(F(2)**-24,F(2)**-20,enlargement=2),dict(first=19,last=25,count=7))
    def test_strict_stage_endpoint_retained(self):self.assertEqual(universal_bounds().stage_bound(F(2)**-24,F(2)**-16),7)
    def test_closed_terminal_mesh_centers(self):self.assertEqual(universal_bounds().closed_centers(0,0,F(1,10**36)),(-1,1))
    def test_same_next_time_curl_requirement(self):
        for m in (0,1,5,100):self.assertEqual(universal_bounds().requirements(m)['potential_for_velocity_Cauchy'],m+2)
    def test_pressure_Cauchy_requirement(self):self.assertEqual(universal_bounds().requirements(7)['pressure_Cauchy'],8)
    def test_no_inverse_eta_endpoint_distance(self):self.assertIn('no 1/(1-eta^2)',audit.protocol()['closed_atlas'])
    def test_complete_axis_seed_curl_components(self):self.assertEqual(universal_checks().exact_review()['exact_seed_curl_components'],3)
    def test_all_36_mixed_pressure_identities(self):self.assertEqual(len(universal_bounds().mixed_review()['exact_terminal_rows']),36);self.assertTrue(all(universal_bounds().mixed_review()['checks'].values()))
    def test_heat_factorial_bounds_at_all_control_orders(self):
        for row in gevrey_review()['rows']:self.assertTrue(row['checks']['heat_lower'] and row['checks']['heat_upper'])
    def test_pressure_factorial_bounds_full_convolution(self):
        for row in gevrey_review()['rows']:self.assertTrue(row['checks']['pressure_lower'] and row['checks']['pressure_upper'])
    def test_reciprocal_binomial_sum_bound(self):self.assertTrue(all(r['checks']['reciprocal_binomial_sum_at_most_three'] for r in gevrey_review()['rows']))
    def test_pressure_divisor_not_discarded_in_sharp_lower_bound(self):self.assertTrue(gevrey_row(F(1,200),64)['checks']['pressure_divisor_retained'])
    def test_first_derivative_exact_endpoint_factors(self):
        b=universal_bounds();h=F(1,200);self.assertEqual(b.heat_A(h,1),h*(1+h));self.assertEqual(b.heat_S(h,1),2*h*(1+h))
    def test_internal_pressure_cross_cannot_be_omitted(self):
        b=universal_bounds();h=F(1,200);self.assertGreater(b.heat_S(h,2),2*b.heat_A(h,2))
    def test_factorial_bound_rejects_zero_derivative_for_lower_formula(self):
        with self.assertRaises(ValueError):gevrey_row(F(1,200),0)
    def test_original_h_domain_not_extended(self):
        for h in (0,-1,F(1,100)):
            with self.assertRaises(ValueError):gevrey_row(h,2)
    def test_Gevrey_not_promoted_to_full_interior(self):self.assertFalse(audit.FLAGS['full_canonical_interior_Gevrey_order_or_analyticity_claimed'])
    def test_nonzero_exterior_qualification_recorded(self):self.assertIn('C!=0',audit.protocol()['Gevrey_scope'])
    def test_uniform_exterior_rate_and_all_order_scope(self):self.assertTrue(gevrey_review()['checks']['pure_time_uniform_rate_is_eight_over_rmin_squared']);self.assertTrue(gevrey_review()['checks']['all_order_proof_not_finite_row_scan'])
    def test_exact_oscillatory_first_time_derivative(self):self.assertEqual(dt({(0,'sin'):F(1)}),{(-2,'cos'):F(1)})
    def test_exact_oscillatory_second_time_derivative(self):self.assertEqual(dt(dt({(0,'sin'):F(1)})),{(-3,'cos'):F(2),(-4,'sin'):F(-1)})
    def test_all_nine_differentiated_exact_residual_cancellations(self):self.assertEqual(len(obstruction_review()['rows']),9);self.assertTrue(all(r['exact_momentum_cancellation'] for r in obstruction_review()['rows']))
    def test_omitting_counterexample_pressure_is_detected(self):self.assertTrue(dt({(0,'sin'):F(1)}))
    def test_invalid_quadrature_rejected(self):
        with self.assertRaises(ValueError):dt({(0,'bad'):F(1)})
    def test_counterexample_not_actual_canonical_field(self):self.assertFalse(obstruction_review()['actual_canonical_field'])
    def test_two_distinct_terminal_subsequences_recorded(self):self.assertIn('+1 and -1',obstruction_review()['unequal_terminal_subsequences'])
    def test_independent_positive_Gamma_heat_comparisons(self):self.assertEqual(universal_checks().integral_review()['nonterminal_heat_comparisons'],24)
    def test_independent_complete_pressure_product_comparisons(self):self.assertEqual(universal_checks().integral_review()['nonterminal_pressure_product_comparisons'],18)
    def test_all_finite_Taylor_and_Cauchy_controls(self):self.assertEqual(universal_checks().exact_review()['finite_Taylor_controls'],21);self.assertEqual(universal_checks().exact_review()['Cauchy_controls'],21)
    def test_frozen_record_exact_roundtrip(self):
        f=json.loads((audit.HERE/'evidence.json').read_text());audit.validate(f);audit.compare_records(self.record,f);audit.compare_records(self.record,json.loads(json.dumps(self.record)))
    def test_forbidden_unconditional_future_and_origin_promotions_rejected(self):
        for flag in ('unconditional_complete_residual_or_blowup_proved','regular_velocity_or_pressure_at_singular_origin_inferred','future_heat_integral_analytic_continuation_or_zero_pressure_gluing_licensed','full_canonical_interior_Gevrey_order_or_analyticity_claimed'):
            bad=deepcopy(self.record);bad['flags'][flag]=True
            with self.assertRaises(ValueError):audit.validate(bad)
    def test_source_hash_tamper_rejected(self):
        bad=deepcopy(self.record);key=next(iter(bad['source_hashes']));bad['source_hashes'][key]='0'*64
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_same_count_provenance_substitution_rejected(self):
        bad=deepcopy(self.record);v=bad['inherited_provenance']['bindings'];v['forged']=v.pop(next(iter(v)))
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_factorial_review_tampering_rejected(self):
        bad=deepcopy(self.record);bad['reviews']['sharp_exterior_Gevrey_two']['rows'][0]['heat_normalized']='0'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_cached_reviews_not_mutated(self):
        bad=deepcopy(self.record);bad['reviews']['exact_residual_endpoint_obstruction']['checks'].clear()
        self.assertTrue(audit.reviews()['exact_residual_endpoint_obstruction']['checks'])
    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'kept.json';p.write_text('keep')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
if __name__=='__main__':unittest.main()
