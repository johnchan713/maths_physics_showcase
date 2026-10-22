#!/usr/bin/env python3
"""Implemented cutoff, all-order recurrences and unchanged-selector scope checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import gaussian_polynomial,step_bounds,cutoff_bounds,verify_schedule,ceil_log2,comparison_value,selector_gate,fixture
from controls import review,supplied_stages


class BorelCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate();cls.f=fixture()

    def control(self,part,fragment):
        rows={k:v for k,v in self.result['parts'][part]['checks'].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_frozen_implemented_cutoff_and_squared_flat_step_are_identified(self):
        self.control('source','implemented_step_body_is_the_squared_flat_logistic')
        self.control('source','implemented_Borel_cutoff_is_exact_same_affine_step')

    def test_floating_cutoff_helper_is_not_a_norm_certificate(self):
        self.control('source','prior_numeric_cutoff_helper_is_not_used_as_an_interval_certificate')

    def test_exact_Gaussian_polynomial_recurrence(self):
        self.control('jets','Gaussian_raw_recurrence_identity')
        self.assertEqual(gaussian_polynomial(1),{3:2})
        self.assertEqual(gaussian_polynomial(2),{4:-6,6:4})

    def test_all_fixed_order_constants_and_raw_reciprocal_rule(self):
        self.control('recurrences','all_fixed_step_bounds_finite_positive')
        self.control('recurrences','raw_binomial_reciprocal_rule_retained')
        self.assertEqual(step_bounds(4)['reciprocal'][0],81)

    def test_low_sharp_bounds_improve_the_independent_all_order_route(self):
        self.control('recurrences','sharp_low_bounds_do_not_exceed_independent_coarse_route')
        self.assertGreater(step_bounds(2)['coarse'][2],1216)

    def test_Euler_norms_are_padded_through_order(self):
        self.control('recurrences','padded_Euler_norms_are_monotone')
        self.assertEqual(cutoff_bounds(2)['K'],[1,64,4928])

    def test_second_Euler_jet_retains_first_ordinary_derivative(self):
        self.control('jets','ordinary_second_jet_cannot_replace_Euler_second')
        self.control('recurrences','Euler_second_requires_first_and_second_ordinary_cutoff_rows')

    def test_actual_implemented_midpoint_rules_out_toy_K2_three(self):
        self.control('jets','midpoint_first_and_second_Euler_jets_are_minus_twelve')
        self.control('jets','padded_K2_three_cannot_bound_this_implemented_cutoff')

    def test_every_checked_midpoint_Euler_jet_is_bounded(self):
        self.control('jets','exact_Euler_midpoint_jet_bounded')

    def test_qualitative_plateau_conditions_alone_do_not_bound_norms(self):
        self.control('jets','qualitative_plateaus_do_not_bound_cutoff_norms')

    def test_supplied_complete_schedule_prefix_passes_without_mutation(self):
        self.control('schedule','all_supplied_normalized_prefix_rows_pass')
        self.control('schedule','no_existing_scale_is_selected_or_changed')

    def test_missing_normalized_row_or_nonconsecutive_index_is_rejected(self):
        stages=supplied_stages();stages[-1]['constants']=stages[-1]['constants'][:-1]
        with self.assertRaises(ValueError):verify_schedule(stages)
        stages=supplied_stages();stages[1]['n']=3
        with self.assertRaises(ValueError):verify_schedule(stages)

    def test_bad_unchanged_scale_reports_failure_instead_of_reselection(self):
        self.control('schedule','inadequate_existing_scale_reports_normalized_row_failure')
        self.control('schedule','insufficient_separation_is_not_silently_repaired')

    def test_existing_q_domain_gate_is_retained(self):
        self.control('schedule','original_common_q_domain_gate_is_retained')

    def test_finite_prefix_is_not_an_infinite_schedule_certificate(self):
        self.control('schedule','prefix_is_not_infinite_schedule_certification')
        self.assertFalse(verify_schedule(supplied_stages())['infinite_schedule_certified'])
        self.assertFalse(cutoff_bounds(6)['universal_retained_schedule_verified'])

    def test_binary_log_upper_is_exact_on_rational_edges(self):
        self.control('schedule','integer_binary_log_upper')
        self.assertEqual(ceil_log2(2**500),500)
        self.assertEqual(ceil_log2(Q(2**500+1)),501)
        self.assertEqual(ceil_log2(Q(2**500-1)),500)

    def test_complete_seven_group_Cphi_matches_piecewise_affine_cutoff_recipe(self):
        self.control('comparison','all_seven_groups_equal_cutoff_dependence')

    def test_zero_first_tuple_norm_keeps_saturated_branch(self):
        r=self.f['radial'];o=r['result'];from budgets import comparison_dependence
        dep=comparison_dependence(r['h'],o['leading']['through'][1],o['leading']['through'][2],0,r['geometry'])
        self.assertEqual(comparison_value(dep,1),comparison_value(dep,4928))

    def test_radial_plus_axial_cutoff_cost_has_factor_five(self):
        dep=self.f['dependence'];r=self.f['radial'];o=r['result'];g=r['geometry'];from budgets import head
        c=head.chart_constants(1,g['Rmin'],g['Rmax']);car=o['comparison']['carriers']
        self.assertEqual(dep['E'],5*g['Lhigh']*(car['theta']+car['axial'])*(c['F']+c['G']))

    def test_certified_K2_joins_same_full_pressure_first_repair(self):
        self.control('comparison','certified_implemented_cutoff_bound_joins_same_radial_repair')
        self.control('comparison','earlier_K_three_is_only_a_manufactured_bound')

    def test_necessary_and_sufficient_recipe_tests_remain_distinct(self):
        self.control('comparison','midpoint_lower_recipe_can_block_all_admissible_cutoff_upper_norms')
        self.control('comparison','large_supplied_length_can_pass_both_recipe_gates')
        out=self.f['gate']
        self.assertLess(out['sufficient_recipe_margin'],out['necessary_recipe_margin'])

    def test_inconsistent_growth_and_tilt_selector_data_are_rejected(self):
        dep=self.f['dependence']
        with self.assertRaises(ValueError):selector_gate(dep,4,2,1,100000)
        with self.assertRaises(ValueError):selector_gate(dep,1,1,1,100000)

    def test_negative_nonexact_or_unpadded_inputs_are_rejected(self):
        with self.assertRaises(ValueError):gaussian_polynomial(-1)
        with self.assertRaises(ValueError):ceil_log2(Q(1,2))
        with self.assertRaises(TypeError):ceil_log2(1.0)
        with self.assertRaises(ValueError):comparison_value(self.f['dependence'],Q(1,2))
        dep=dict(self.f['dependence']);dep['D']=-1
        with self.assertRaises(ValueError):comparison_value(dep,64)
        stages=supplied_stages();stages[0]['gain_lower']=0
        with self.assertRaises(ValueError):verify_schedule(stages)

    def test_six_frozen_records_and_all_declared_source_hashes_match(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),6)
        self.assertTrue(all(self.record['gates'].values()))

    def test_new_frozen_evidence_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)

    def test_actual_sequence_selector_and_residual_scope_cannot_be_promoted(self):
        for key,value in self.record['flags'].items():
            if not value:
                changed=deepcopy(self.record);changed['flags'][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)

    def test_remaining_identity_schedule_and_domain_obligations_cannot_be_removed(self):
        changed=deepcopy(self.record);changed['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(changed)
        changed=deepcopy(self.record);changed['source_hashes'].pop(next(iter(changed['source_hashes'])))
        with self.assertRaises(ValueError):audit.validate(changed)

    def test_frozen_output_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
