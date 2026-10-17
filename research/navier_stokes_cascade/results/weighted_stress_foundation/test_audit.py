#!/usr/bin/env python3
"""Weighted edge, kernel endpoint, beta term, exact support and scope controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from edge_bounds import binomial,volterra_coefficients,flat_power_polynomial,cutoff_integer_log,actual_input_review,outer_review,background_review
from edge_checks import maps,volterra_jet,divide_delta,kernel_stats,outer_jets,outer_stats,heat_stats,budget_stats,diagnostic_review
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,main


class WeightedStressFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_actual_leading_parent_and_frozen_parameters(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'a779312a3036c03ca7abea540b18e19785c9bc03')
        self.assertIn('t1=C^-20',self.record['protocol']['unchanged_parameters'])

    def test_retained_exact_core_and_finite_lemmas_not_erased(self):
        self.assertTrue(FLAGS['exact_pressure_family_core_and_finite_existence_inverse_gap_lemmas_remain_dependencies'])
        self.assertIn('identically zero',actual_input_review()['exact_core_requirement'])

    def test_actual_activated_not_unmodified_core_used_at_inner_edge(self):
        self.assertTrue(actual_input_review()['checks']['actual_activation_not_extended_unmodified_core_used'])
        self.assertIn('(1-kappa0)',actual_input_review()['inner_factor'])

    def test_flat_Volterra_endpoint_constant_and_first_nontrivial_Jacobian_jet(self):
        self.assertEqual(volterra_coefficients([3,5,7,11],4),[3,5,7,11-Q(5,8)])
        self.assertEqual(volterra_coefficients([1,0,0,0],4),[1,0,0,0])

    def test_all_768_kernel_and_independent_ODE_endpoint_components(self):
        r=kernel_stats();self.assertEqual(r['count'],768);self.assertTrue(r['match'])

    def test_fractional_binomial_radial_Jacobian_not_replaced_by_one(self):
        self.assertEqual(binomial(-Q(1,2),2),Q(3,8))
        self.assertEqual(kernel_stats()['missing_jacobian_coefficient'],-Q(1,8))

    def test_mixed_parameter_kernel_jets_match_independent_ODE(self):
        self.assertTrue(kernel_stats()['mixed_coefficients_equal'])
        m=maps();j=m.jet({(1,2):Q(3),(3,1):Q(5)},2,8)
        self.assertFalse((volterra_jet(j)-j).terms=={})

    def test_kernel_rejects_wrong_constant_empty_jets_or_method(self):
        for c,j,method in ((0,[1],'kernel'),(-1,[1],'ode'),(4,[],'ode'),(4,[1],'Taylor')):
            with self.assertRaises(ValueError):volterra_coefficients(j,c,method)

    def test_delta_division_requires_exact_vanishing_not_small_numbers(self):
        m=maps();s=m.variable(0,order=8)
        self.assertFalse((divide_delta(s**3*(1+s),3)-(1+s)).terms)
        with self.assertRaises(ValueError):divide_delta(s**3+Q(1,10**100),3)

    def test_flat_derivatives_keep_inner_constant_and_outer_inverse_power(self):
        self.assertEqual(flat_power_polynomial(Q(1,10000),0,1),{3:Q(1,5000)})
        self.assertEqual(flat_power_polynomial(4,-3,1),{3:Q(8),1:Q(-3)})
        self.assertNotEqual(flat_power_polynomial(1,-3,1),flat_power_polynomial(4,-3,1))

    def test_second_flat_derivative_retains_all_logarithmic_derivative_terms(self):
        self.assertEqual(flat_power_polynomial(4,0,2),{6:Q(64),4:Q(-24)})
        for n in range(16):self.assertEqual(max(flat_power_polynomial(4,-3,n)),3*n)

    def test_all_675_factored_outer_mixed_components_match_kernel_ODE(self):
        r=outer_stats();self.assertEqual(r['count'],675);self.assertTrue(r['equal'])

    def test_actual_finite_Gamma_pole_jets_have_positive_angular_factor(self):
        for e in (Q(-1),Q(1)):
            a=outer_jets(e);self.assertEqual(a['b0'].value(),1)
            self.assertNotEqual(a['bz'].value(),0)
        self.assertFalse(FLAGS['heat_pole_Taylor_series_assumed_convergent'])

    def test_outer_endpoint_quotient_keeps_large_Xb_factor(self):
        r=outer_stats();self.assertTrue(r['endpoints'])
        self.assertEqual(outer_jets(Q(1))['W'].value(),Q(256)/Q(99,100)**2)
        self.assertEqual(outer_jets(Q(-1))['W'].value(),outer_jets(Q(1))['W'].value())

    def test_eta_zero_first_two_radial_jets_vanish_but_third_does_not(self):
        r=outer_stats();self.assertEqual((r['zero_first'],r['zero_second']),(0,0))
        self.assertEqual(r['cubic'],Q(3,4)*128)

    def test_eta_zero_second_parameter_jet_is_four_Xb(self):
        self.assertEqual(outer_stats()['eta_second'],4*128)

    def test_missing_beta_term_erases_required_cubic(self):
        self.assertEqual(outer_stats()['omitted_beta_cubic'],0)
        self.assertNotEqual(outer_stats()['cubic'],0)

    def test_actual_eta_zero_Gamma_cubic_limit_refines(self):
        r=heat_stats();self.assertTrue(all(r['checks'].values()));self.assertEqual(len(r['rows']),4)

    def test_actual_coefficient_inputs_now_use_closed_leading_foundation(self):
        self.assertTrue(FLAGS['actual_leading_inputs_to_formal_coefficient_maps_supplied'])
        self.assertIn('analytic FIRST collar',actual_input_review()['leading_input'])

    def test_common_higher_stress_support_does_not_include_outer_edge(self):
        s=actual_input_review()['support'];self.assertEqual(s['stress_right'],{'1':'Xb','2':'Xv','n>=3':'Xvel'})
        self.assertIn('tc/16',s['inner'])

    def test_first_two_pressure_tails_and_leading_radial_flow_not_dropped(self):
        self.assertTrue(FLAGS['full_pressure_memory_I4_first_two_orders_and_radial_velocity_retained'])
        self.assertEqual(actual_input_review()['support']['pressure_right']['2'],'Xv')
        self.assertIn('O_I(q^h)',background_review()['radial_error'])

    def test_conservative_backward_identification_keeps_shifted_viscosity(self):
        r=actual_input_review()['checks']
        self.assertTrue(r['exact_five_moments_and_shifted_viscosity_identify_signed_backward_stresses'])
        self.assertTrue(r['order_one_angular_viscosity_uses_renormalized_leading_swirl_integral'])

    def test_one_schedule_includes_all_seminorms_through_each_n(self):
        r=budget_stats();self.assertTrue(r['all_smaller_q']);self.assertEqual(len(r['rows']),13)
        self.assertEqual(r['count'],4*sum(n+1 for n in range(2,15)))

    def test_exact_integer_log_budget_controls_arbitrarily_smaller_q(self):
        h=Q(1,200);ell=cutoff_integer_log(h,2,[10,20,30],previous=50,boundary=10)
        for extra in (0,1,1000000):self.assertGreater(2*h*(ell+extra),32)

    def test_budget_rejects_omitted_derivative_order_or_h_zero(self):
        for h,n,logs in ((0,2,[1,2,3]),(Q(1,200),1,[1,2]),(Q(1,200),2,[1,2])):
            with self.assertRaises(ValueError):cutoff_integer_log(h,n,logs)

    def test_generic_weighted_smallness_cannot_replace_relative_quotient(self):
        self.assertLess(1-Q(1,100)/Q(1,200),0)
        self.assertTrue(FLAGS['T1_over_abs_T0_smooth_on_closed_annulus_including_eta_zero'])

    def test_direct_stress_sum_is_not_zero_complete_NS_residual(self):
        self.assertTrue(background_review()['checks']['direct_stress_sum_not_identified_with_zero_complete_NS_residual'])
        self.assertFalse(FLAGS['unconditional_complete_residual_improvement_verified'])

    def test_frozen_record_and_JSON_roundtrip_match_sources(self):
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.record,f)
        validate(json.loads(json.dumps(self.record)))

    def test_all_836_inherited_and_19_direct_source_bindings_preserved(self):
        self.assertEqual(len(self.record['inherited_provenance']),836)
        self.assertEqual(len(self.record['source_hashes']),19)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_cached_diagnostic_reviews_are_snapshots(self):
        r=generate();key='omitting_beta_erases_required_nonzero_cubic'
        r['reviews']['independent_exact_kernel_heat_and_budget_controls']['checks'][key]=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(diagnostic_review()['checks'][key])

    def test_provenance_key_substitution_rejected_even_with_same_count(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])

    def test_analyticity_uniformity_and_unconditional_scope_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)

    def test_endpoint_or_pressure_support_formula_tampering_rejected(self):
        for key in ('eta_zero_jets','relative_quotient'):
            r=deepcopy(self.record);r['reviews']['all_order_outer_factor_and_eta_zero_cubic'][key]='0'
            with self.assertRaises(ValueError):validate(r)
        r=deepcopy(self.record);r['reviews']['actual_inner_factor_and_coefficient_support_inputs']['support']['pressure_right']['2']='Xvel'
        with self.assertRaises(ValueError):validate(r)

    def test_remaining_gates_bindings_and_actual_pulse_obligations_not_erased(self):
        for key in ('remaining_obligations','gates','source_hashes'):
            r=deepcopy(self.record)
            if key=='remaining_obligations':r[key].pop()
            else:r[key].pop(next(iter(r[key])))
            with self.assertRaises(ValueError):validate(r)

    def test_cli_refuses_overwriting_existing_evidence_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')


if __name__=='__main__':unittest.main()
