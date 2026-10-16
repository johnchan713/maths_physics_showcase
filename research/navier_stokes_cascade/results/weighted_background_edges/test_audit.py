#!/usr/bin/env python3
"""Weighted edge, relative quotient, cutoff quantifier and scope regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from audit import HERE,FLAGS,compare_records,generate,main,validate
from edges import flat_polynomial,global_log_jet,log_flat_majorant,stress_loss,tail_split


class WeightedEdgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=generate()

    def test_all_exact_flat_polynomial_controls(self):
        r=self.record['reviews']['exact_flat_algebra_and_support_ledger']
        self.assertEqual(r['exact_polynomial_checks'],75)
        self.assertTrue(all(r['checks'].values()))

    def test_zero_and_first_flat_derivatives(self):
        self.assertEqual(flat_polynomial(Q(1,10000),0),{0:Q(1)})
        self.assertEqual(flat_polynomial(Q(1,10000),1),{3:Q(1,5000)})

    def test_second_derivative_keeps_chain_rule_term(self):
        self.assertEqual(flat_polynomial(4,2),{6:Q(64),4:Q(-24)})

    def test_high_order_degree_is_three_times_radial_order(self):
        self.assertEqual(max(flat_polynomial(Q(1,10**6),25)),75)

    def test_invalid_orders_and_exponential_constants_rejected(self):
        for c,n in ((0,1),(-1,1),(1,-1),(1,True),(1,1.5)):
            with self.assertRaises(ValueError):
                flat_polynomial(c,n)

    def test_both_edge_derivative_signs_are_retained(self):
        self.assertEqual(global_log_jet(Q(1),Q(4),Q(1,2),Q(1,4),1),Q(16)-Q(512))

    def test_total_order_weight_loss_includes_outer_pole(self):
        self.assertEqual(stress_loss(0),3)
        self.assertEqual(stress_loss(10),33)

    def test_tail_starts_late_enough_even_for_zero_derivatives(self):
        self.assertEqual(tail_split(0),dict(split=2,tail_gain_multiple=3,coefficient_multiple=Q(1,4),retained_gain_multiple=2))

    def test_arbitrary_fixed_derivative_uses_same_two_h_gain(self):
        self.assertEqual(tail_split(100)['retained_gain_multiple'],2)
        self.assertEqual(tail_split(100)['tail_gain_multiple'],101)

    def test_two_weight_constants_match_independent_differentiation(self):
        r=self.record['reviews']['both_edge_weights_and_fixed_derivatives']
        self.assertLess(r['maximum_flat_derivative_error'],1e-65)
        self.assertLess(r['maximum_two_edge_error'],1e-65)

    def test_exponential_dominates_any_fixed_inverse_power(self):
        r=self.record['reviews']['both_edge_weights_and_fixed_derivatives']
        self.assertEqual(r['flat_majorant_checks'],180)
        self.assertTrue(r['checks']['exponential_absorbs_every_fixed_polynomial_loss'])

    def test_majorant_critical_point_is_exact_maximum(self):
        with mp.workdps(50):
            c=mp.mpf(4);r=mp.mpf(12);y=mp.sqrt(2*c/r)
            self.assertLess(abs(-c/y**2-r*mp.log(y)-log_flat_majorant(c,3,9)),mp.mpf('1e-45'))

    def test_invalid_majorants_rejected(self):
        for c,N,p in ((0,1,0),(-1,1,0),(1,-1,0),(1,1,-1)):
            with self.assertRaises(ValueError):
                log_flat_majorant(c,N,p)

    def test_replacing_width_squared_by_one_is_detected(self):
        r=self.record['reviews']['both_edge_weights_and_fixed_derivatives']
        self.assertGreater(r['wrong_weight_log_excess'],1e9)
        self.assertEqual(self.record['protocol']['inner_exponent'],'t1^2/ya^2')

    def test_slower_weight_also_breaks_positive_lower_bound(self):
        self.assertTrue(self.record['reviews']['both_edge_weights_and_fixed_derivatives']['checks']['slower_weight_defeats_positive_lower_comparison'])

    def test_all_smaller_q_schedule_requirements_pass(self):
        r=self.record['reviews']['all_smaller_q_diagonal_schedule']
        self.assertEqual(len(r['rows']),11)
        self.assertGreater(r['minimum_log_margin'],0)
        self.assertTrue(all(r['checks'].values()))

    def test_both_weighted_and_relative_seminorms_included(self):
        self.assertEqual(self.record['reviews']['all_smaller_q_diagonal_schedule']['additional_seminorms'],['Tn/zeta through n','Tn/abs(T0) through n'])

    def test_cutoff_scales_not_claimed_numerically_evaluated(self):
        self.assertFalse(self.record['protocol']['actual_cutoff_constants_or_scales_numerically_evaluated'])

    def test_factored_outer_quotient_has_no_raw_flat_division(self):
        r=self.record['reviews']['smooth_relative_heat_edge_and_direction']
        self.assertLess(r['maximum_factored_identity_error'],1e-55)
        self.assertTrue(all(r['checks'].values()))

    def test_parameter_endpoints_and_zero_eta_are_present(self):
        rows=self.record['reviews']['smooth_relative_heat_edge_and_direction']['rows']
        self.assertEqual({row['eta'] for row in rows},{'-1.0','0.0','0.5','1.0'})
        self.assertEqual(len(rows),12)

    def test_relative_smallness_bounds_norm_and_direction(self):
        for row in self.record['reviews']['smooth_relative_heat_edge_and_direction']['rows']:
            self.assertLessEqual(row['norm_ratio'],1.5)
            self.assertGreaterEqual(row['norm_ratio'],.5)
            self.assertLessEqual(row['new_cone_ratio'],row['leading_cone_ratio']+1e-60)

    def test_generic_weighted_bound_does_not_preserve_cone(self):
        self.assertEqual(self.record['reviews']['exact_flat_algebra_and_support_ledger']['weighted_counterexample_total_sign'],'-1')

    def test_same_fixed_N_for_every_finite_jet(self):
        self.assertTrue(self.record['protocol']['finite_N_not_reselected_for_higher_derivatives'])
        self.assertTrue(self.record['protocol']['common_q_domain_for_every_fixed_jet'])

    def test_pressure_tails_not_removed_by_new_stress_estimate(self):
        self.assertTrue(FLAGS['first_two_reserved_mean_pressure_tails_retained'])

    def test_parent_and_all_564_historical_hashes_replayed(self):
        self.assertTrue(self.record['gates']['exact_parent_supports_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),564)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_record_reproduces(self):
        frozen=json.loads((HERE/'evidence.json').read_text())
        validate(frozen);compare_records(self.record,frozen)

    def test_pulse_mean_residual_blowup_and_numerical_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):
                    validate(r)

    def test_analytic_lemma_dependency_cannot_be_removed(self):
        r=deepcopy(self.record);r['flags']['relies_on_replayed_core_loop_repair_support_and_heat_lemmas']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_shortening_edge_weight_loss_or_tail_split_rejected(self):
        for field,value in (('total_order_weight_loss','0'),('inner_exponent','1/ya^2')):
            r=deepcopy(self.record);r['protocol'][field]=value
            with self.assertRaises(ValueError):
                validate(r)

    def test_missing_gate_or_obligation_rejected(self):
        for section in ('gates','remaining_obligations'):
            r=deepcopy(self.record)
            if section=='gates':
                r[section].pop(next(iter(r[section])))
            else:
                r[section].pop()
            with self.assertRaises(ValueError):
                validate(r)

    def test_changed_provenance_rejected(self):
        r=deepcopy(self.record);r['inherited_provenance'][next(iter(r['inherited_provenance']))]=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_changed_current_source_rejected(self):
        r=deepcopy(self.record);r['source_hashes'][next(iter(r['source_hashes']))]='0'*64
        with self.assertRaises(ValueError):
            validate(r)

    def test_nonfinite_diagnostic_rejected(self):
        r=deepcopy(self.record);r['reviews']['all_smaller_q_diagonal_schedule']['minimum_log_margin']=float('nan')
        with self.assertRaises(ValueError):
            compare_records(self.record,r)

    def test_CLI_cannot_overwrite_frozen_record(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':
    unittest.main()
