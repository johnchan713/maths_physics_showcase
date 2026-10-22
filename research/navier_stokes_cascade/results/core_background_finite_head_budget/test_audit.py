#!/usr/bin/env python3
"""Finite-head bounds, exact manufactured jets, original-selector and scope guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import BackgroundBudget,P,phase,fixture,finite_head,increment_constant,chart_constants,length_compatibility
from controls import review


class BackgroundHeadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def control(self,part,fragment):
        rows={k:v for k,v in self.result['parts'][part]['checks'].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_finite_head_keeps_first_coefficient_at_low_orders(self):
        self.control('tail','coefficient_one_always_in_finite_head')
        self.assertEqual(finite_head(0,{1:Q(7)})['correction'],Q(15,2))
        with self.assertRaises(ValueError):finite_head(0,{})

    def test_all_q_geometric_tail_has_exact_closed_sum(self):
        self.control('tail','entire_geometric_tail_bound')
        self.control('tail','whole_interval_correction_constant')

    def test_weak_diagonal_from_one_cannot_give_two_h_rate(self):
        self.control('tail','starting_diagonal_tail_at_one_loses_q_to_2h_rate')

    def test_head_inputs_end_at_the_correct_order(self):
        trace=set();b=fixture(trace);b.bound(5)
        self.assertEqual({t[1] for t in trace if t[0]=='profile'},{1,2,3,4})
        self.assertEqual({t[2] for t in trace if t[0]=='profile'},{5})
        self.assertIn(('cutoff',6),trace)

    def test_radial_cutoff_curl_commutator_not_deleted_or_reversed(self):
        self.control('cutoff','radial_cutoff_curl_commutator_nonzero')
        self.control('cutoff','commutator_sign_is_essential')

    def test_Euler_product_rule_and_extra_cutoff_jet(self):
        self.control('cutoff','Euler_product_rule_retains_commutator')
        self.control('cutoff','cutoff_extra_Euler_jet_bound')
        self.assertEqual(fixture().bound(3)['cutoff_input_order'],4)

    def test_frozen_q_derivatives_cannot_be_substituted_unchanged(self):
        self.control('cutoff','Euler_second_jet_differs_from_frozen_second_jet')
        with self.assertRaises(ValueError):BackgroundBudget(Q(1,200),None,None,None,1,2,'frozen')

    def test_chart_retains_all_three_normalized_velocity_components(self):
        self.control('chart','F_finite_chart_composition_bounded')
        self.control('chart','G_finite_chart_composition_bounded')
        self.control('chart','btilde_finite_chart_composition_bounded')

    def test_chart_uses_positive_annular_reciprocal(self):
        self.control('chart','radial_reciprocal_is_not_dropped')
        self.assertEqual(chart_constants(0,1,2)['F'],4)
        self.assertEqual(chart_constants(0,1,2)['btilde'],2)
        with self.assertRaises(ValueError):chart_constants(2,0,2)

    def test_complete_correction_enters_all_background_bounds(self):
        b=fixture();larger=BackgroundBudget(b.h,b.leading,lambda n,m:10*b.profiles(n,m),b.cutoff,1,2,'euler')
        for name in ('F','G','btilde','error_FG','full_FG'):
            self.assertGreater(larger.bound(1)[name],b.bound(1)[name])
        self.assertEqual(larger.leading_bound(2),b.leading_bound(2))

    def test_error_normalization_keeps_s_to_two_h_factor(self):
        b=fixture();out=b.bound(2)
        self.assertEqual(out['error_FG'],4*(out['chart']['F']+out['chart']['G'])*out['head']['correction'])

    def test_frozen_shear_cancellation_is_exact(self):
        self.control('comparison','frozen_shear_cancellation')

    def test_normal_errors_retain_full_path_rounding_and_axial_terms(self):
        self.control('comparison','full_normal_error_envelope')
        groups=fixture().comparison(1,1,1,1,2,1)['groups']
        self.assertGreater(groups['radial_rounding'],0)
        self.assertGreater(groups['axial_transport'],0)

    def test_seven_group_comparison_uses_complete_positive_cost(self):
        self.control('comparison','Cphi_sums_all_seven_groups')
        self.control('comparison','every_error_group_has_positive_cost')
        out=fixture().comparison(1,1,1,1,2,1)
        for value in out['groups'].values():self.assertLess(out['Cphi']-value,out['Cphi'])

    def test_comparison_requires_only_finite_stated_lower_jet_set(self):
        self.control('comparison','comparison_uses_only_coefficient_one_jets_through_one')
        self.control('comparison','comparison_uses_leading_two_and_cutoff_two')

    def test_existing_selector_checks_remain_in_force(self):
        self.control('selector','unchanged_selector_all_original_checks_pass')
        self.control('selector','slow_gate_already_dominates_parent_rounding_gate')

    def test_original_length_margin_joins_parent_value_gate(self):
        self.control('selector','unchanged_length_budget_joins_parent_normal_margin')
        self.control('selector','normal_floor_same_as_parent')
        out=length_compatibility(1,2,1,Q(1,10),12,1,1)
        self.assertEqual(out['beta'],Q(1,100));self.assertEqual(out['inverse_normal_square'],400)

    def test_failed_length_budget_is_not_silently_reselected(self):
        self.control('selector','insufficient_original_length_is_reported_not_reselected')
        self.assertFalse(length_compatibility(1,2,100,1,1,1,1)['gates']['original_length_budget_dominates_normal_margin'])

    def test_explicit_background_values_feed_actual_phase_source_coefficients(self):
        b=fixture();p=b.phase_inputs(lambda n:P({1:2}) if n==0 else P(int(n==1)),1,2,3,Q(1,2))
        values=p.source_coefficients(2)
        self.assertEqual(set(values),{'b','G','V','r','normal','K','Eik','bR','Bz'})
        self.assertTrue(all(v.evaluate(1)>0 for v in values.values()))
        self.assertEqual(p.data('btilde',2).evaluate(1),b.bound(2)['btilde'])

    def test_explicit_jets_feed_full_pressure_inverse_with_remaining_premises(self):
        p=fixture().phase_inputs(lambda n:P({1:2}) if n==0 else P(int(n==1)),1,2,3,Q(1,2))
        gate=phase.normal_gate(1,2,1,100,Q(1,20000))
        forced=p.forced(lambda j:P(1),P(2),P({1:2}),gate)
        self.assertGreater(forced.velocity(0).evaluate(1),0)
        self.assertGreater(forced.pressure(1).evaluate(1),0)
        self.assertEqual(forced.profile_rate(3).c,p.profile_rate(3).c)

    def test_bounds_are_finite_exact_rationals_with_no_float_inputs(self):
        for m in range(5):
            for k in ('F','G','btilde'):self.assertIsInstance(fixture().bound(m)[k],Q)
        with self.assertRaises(TypeError):increment_constant(.005,1,1,1,1)
        with self.assertRaises(TypeError):chart_constants(1,1.0,2)

    def test_invalid_profile_cutoff_order_and_geometry_are_rejected(self):
        for args in ((Q(1,200),0,1,1,1),(Q(1,100),1,1,1,1),(Q(1,200),1,1,-1,1),(Q(1,200),1,1,1,0)):
            with self.assertRaises(ValueError):increment_constant(*args)
        with self.assertRaises(ValueError):fixture().bound(-1)
        with self.assertRaises(ValueError):fixture().comparison(1,1,0,1,2,1)
        with self.assertRaises(ValueError):length_compatibility(1,2,1,1,1,Q(1,2),1)

    def test_six_historical_records_and_every_declared_source_match(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),6)
        self.assertTrue(all(self.record['inherited_source_checks'].values()))
        self.assertTrue(all(self.record['gates'].values()))

    def test_new_evidence_reproduces_exactly(self):
        self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)

    def test_actual_constants_selector_and_residual_cannot_be_promoted(self):
        for key,value in self.record['flags'].items():
            if not value:
                changed=deepcopy(self.record);changed['flags'][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)

    def test_remaining_obligation_and_source_binding_cannot_be_deleted(self):
        for field in ('source_hashes','remaining_obligations'):
            changed=deepcopy(self.record)
            if isinstance(changed[field],dict):changed[field].pop(next(iter(changed[field])))
            else:changed[field].pop(0)
            with self.assertRaises(ValueError):audit.validate(changed)

    def test_existing_frozen_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
