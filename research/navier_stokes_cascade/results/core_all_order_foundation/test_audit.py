#!/usr/bin/env python3
"""All-order local core bounds, exact coordinate jets and scope controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from math import comb,factorial
from pathlib import Path
import tempfile
import unittest
from foundation import mixed_bound,analytic_review,zero_free_review,derivative_weight,physical_review,profile_review
from coordinate_checks import Jet,physical_field,normalized_field,operator_word,coordinate_stats
from audit import HERE,FLAGS,generate,validate,compare_records,main,core_review,parent_review


class CoreFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_infinite_core_and_current_parent_reproduced(self):
        self.assertTrue(all(core_review()['checks'].values()))
        self.assertTrue(self.record['gates']['current_conditional_chain_reproduced_without_promotion'])

    def test_joint_geometric_majorant_not_separate_radii(self):
        self.assertEqual(mixed_bound(0,0),Q(8,5))
        with self.assertRaises(ValueError):mixed_bound(0,0,Y=12,eta_offset=Q(3,5))

    def test_mixed_derivative_majorant_against_independent_rational_jet(self):
        Y=Jet.variable(0,5,n=2);eta=Jet.variable(1,Q(1,8),n=2)
        majorant=(1-Y/20-eta)**-1
        for r in range(5):
            for k in range(5-r):self.assertEqual(majorant.partial((r,k)),mixed_bound(r,k))

    def test_mixed_factorial_cannot_be_product_factorials(self):
        r,k=2,3
        wrong=Q(factorial(r)*factorial(k),20**r)/Q(5,8)**(r+k+1)
        self.assertEqual(mixed_bound(r,k)/wrong,comb(r+k,r));self.assertGreater(mixed_bound(r,k),wrong)

    def test_arbitrary_large_fixed_orders_exact_and_finite(self):
        for r,k in ((0,100),(75,40),(100,0)):self.assertGreater(mixed_bound(r,k),0)

    def test_rho_loss_has_actual_parameter_order(self):
        self.assertEqual(mixed_bound(3,7,rho=Q(1,3))/mixed_bound(3,7),3**7)

    def test_invalid_norm_orders_and_margins_rejected(self):
        for r,k in ((-1,0),(True,0),(0,Q(1,2))):
            with self.assertRaises(ValueError):mixed_bound(r,k)
        for args in (dict(rho=0),dict(norm=-1),dict(Y=-1),dict(eta_offset=-1),dict(Y=20)):
            with self.assertRaises(ValueError):mixed_bound(0,0,**args)

    def test_real_bound_and_complex_variation_have_explicit_slack(self):
        r=zero_free_review();self.assertEqual(Q(r['variation_upper']),Q(84,3125))
        self.assertEqual(Q(r['complex_real_part_lower']),Q(741,3125));self.assertGreater(Q(r['complex_real_part_lower']),Q(1,5))

    def test_zero_free_discs_share_one_radius_for_all_orders(self):
        r=zero_free_review();self.assertTrue(all(r['checks'].values()))
        self.assertEqual(r['eta_radius'],'rho/(100*S)');self.assertEqual(r['Y_radius'],'1/(100*S)')

    def test_real_positivity_alone_does_not_give_uniform_complex_strip(self):
        for n in (1,10,10000):
            self.assertEqual(1-n*n*Q(1,n)**2,0)
        self.assertLess(Q(1,10000),Q(1,10))

    def test_pressure_datum_and_increment_budgets_both_present(self):
        r=profile_review();self.assertIn('10*P^2+',r['norm_budgets']['Pi'])
        self.assertIn('p_Y=g^2*Phi^2',r['definitions'])

    def test_radial_average_and_extra_eta_derivative_not_discarded(self):
        r=profile_review();self.assertIn('partial_eta avgU',r['regular_radial_quotient'])
        self.assertTrue(all(r['checks'].values()))

    def test_vanishing_axis_E_is_not_reciprocal_regular(self):
        zero=Jet(0)
        with self.assertRaises(ValueError):zero**-1
        self.assertTrue(profile_review()['checks']['E_equals_radius_times_F_and_cannot_be_inverted_at_axis'])

    def test_exact_implicit_q_jet_contains_high_axial_and_time_orders(self):
        u,p,q=physical_field();self.assertIn((0,0,4,0),q.terms);self.assertIn((0,0,0,4),q.terms)
        self.assertEqual(q.value(),1)

    def test_fractional_jet_power_retains_all_mixed_terms(self):
        x=Jet.variable(0,0,n=2);y=Jet.variable(1,0,n=2)
        f=(1+x+y)**Q(1,2)
        self.assertEqual(f.partial((2,2)),Q(-15,16))
        self.assertEqual(f.partial((1,1)),Q(-1,4))

    def test_exact_all_280_coordinate_identities(self):
        r=coordinate_stats();self.assertEqual(r['multiindices'],70);self.assertEqual(r['identities'],280)
        self.assertEqual(r['maximum_exact_error'],0)

    def test_mixed_coordinate_derivative_orders_commute_with_updated_beta(self):
        r=normalized_field()[3]
        self.assertEqual(operator_word(r,(0,2,3)),operator_word(r,(3,2,0)))

    def test_freezing_beta_breaks_second_time_derivative(self):
        self.assertGreater(coordinate_stats()['frozen_exponent_gap'],0)

    def test_omitting_eta_time_chain_rule_fails(self):
        self.assertGreater(coordinate_stats()['missing_time_eta_gap'],0)

    def test_axis_and_off_axis_divergence_jets_are_zero(self):
        r=coordinate_stats();self.assertEqual(r['axis_divergence'],0);self.assertEqual(r['off_axis_divergence_error'],0)

    def test_axis_radial_quotient_has_nonzero_prescribed_value(self):
        self.assertEqual(coordinate_stats()['axis_radial_derivative'],Q(-699,1798))

    def test_full_pressure_gradient_cannot_be_discarded(self):
        self.assertGreater(coordinate_stats()['pressure_omission_gap'],0)

    def test_full_axial_viscosity_cannot_be_discarded(self):
        self.assertGreater(coordinate_stats()['axial_viscosity_omission_gap'],0)

    def test_manufactured_field_is_not_zero_force_solution(self):
        self.assertGreater(coordinate_stats()['manufactured_full_residual'],0)

    def test_jet_derivative_budget_is_not_silently_exceeded(self):
        z=Jet(1,n=2,order=0)
        with self.assertRaises(ValueError):z.derivative(0)
        with self.assertRaises(ValueError):z.partial((1,0))

    def test_anisotropic_weight_includes_full_time_cost(self):
        self.assertEqual(derivative_weight(2,3,4),1+3*Q(99,200)+4)
        self.assertLess(derivative_weight(0,1,0),derivative_weight(1,0,0))

    def test_invalid_physical_derivative_orders_and_h_rejected(self):
        for args in ((-1,0,0),(1,False,0),(0,0,Q(1,2))):
            with self.assertRaises(ValueError):derivative_weight(*args)
        for h in (0,Q(1,2),-1):
            with self.assertRaises(ValueError):derivative_weight(0,0,0,h)

    def test_actual_unchanged_core_stops_at_four_not_seed_four_point_one(self):
        self.assertEqual(self.record['protocol']['actual_unchanged_leading_interval'],'0<=Y<=4')
        self.assertEqual(self.record['protocol']['seed_real_interval'],'0<=Y<=4.1')
        self.assertIn('<=4;',physical_review()['domain'])

    def test_closed_eta_endpoints_and_positive_q_not_singular_event(self):
        r=physical_review();self.assertTrue(all(r['checks'].values()))
        self.assertIn('singular event excluded',r['domain'])

    def test_758_inherited_hashes_and_13_current_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),758)
        self.assertTrue(all(self.record['inherited_provenance'].values()));self.assertEqual(len(self.record['source_hashes']),13)

    def test_frozen_record_reproduces_and_roundtrips(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r)
        validate(json.loads(json.dumps(self.record)))

    def test_every_unconditional_promotion_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError):validate(r)

    def test_generated_record_does_not_alias_cached_provenance_or_core_controls(self):
        r=generate();c=r['reviews']['infinite_core_premises_reproduced']['checks'];c[next(iter(c))]=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(all(core_review()['checks'].values()))
        r=generate();d=r['inherited_provenance'];d.pop(next(iter(d)));d['invented:binding']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('invented:binding',parent_review()[1])

    def test_precise_majorant_or_pressure_formula_mutation_rejected_behind_gates(self):
        for section,key in (('all_index_joint_analytic_majorant','majorant'),('complete_core_profile_and_axis_quotient','definitions')):
            r=deepcopy(self.record);r['reviews'][section][key]='omit the mixed terms'
            with self.assertRaises(ValueError):validate(r)

    def test_next_activated_collar_obligation_cannot_be_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):validate(r)

    def test_missing_gate_or_source_binding_rejected(self):
        for section in ('gates','source_hashes'):
            r=deepcopy(self.record);r[section].pop(next(iter(r[section])))
            with self.assertRaises(ValueError):validate(r)

    def test_CLI_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':unittest.main()
