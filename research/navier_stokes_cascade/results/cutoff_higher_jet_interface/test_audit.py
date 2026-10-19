"""Exact finite-jet identities, analytical-bound controls and no proof promotion."""
from copy import deepcopy
from fractions import Fraction as F
from math import factorial
import json,tempfile,unittest
from pathlib import Path
import audit
import jet_bounds as b

class HigherCutoffTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_order_zero_identity(self):self.assertEqual(b.log_coefficients(0),(1,));self.assertEqual(b.physical_coefficients(0),{0:1})
    def test_signed_stirling_second_and_third_jacobians(self):self.assertEqual(b.log_coefficients(2),(0,-1,1));self.assertEqual(b.log_coefficients(3),(0,2,-3,1))
    def test_signed_stirling_recurrence_through_ten(self):
        for m in range(1,11):
            v=b.log_coefficients(m);w=b.log_coefficients(m+1)
            for j in range(m+2):self.assertEqual(w[j],(v[j-1] if j else 0)-m*(v[j] if j<len(v) else 0))
    def test_inverse_radius_F_bump_polynomials(self):self.assertEqual(b.log_coefficients(2,F(1,2)),(F(3,4),-2,1))
    def test_all_eighteen_independent_logarithmic_identities(self):self.assertEqual(len(self.r['symbolic']['log_coordinate_comparisons']),18);self.assertTrue(all(v['matches'] for v in self.r['symbolic']['log_coordinate_comparisons']))
    def test_all_nine_physical_chain_identities(self):
        for m in range(9):self.assertTrue(self.r['symbolic']['checks']['physical_chain_order_'+str(m)])
    def test_physical_coefficients_third_and_fourth(self):self.assertEqual(b.physical_coefficients(3),{2:3,3:1});self.assertEqual(b.physical_coefficients(4),{2:3,3:6,4:1})
    def test_physical_coefficients_differentiate_previous_order(self):
        for m in range(1,9):
            old=b.physical_coefficients(m);new=b.physical_coefficients(m+1)
            for k,c in new.items():self.assertEqual(c,(2*k-m)*old.get(k,0)+old.get(k-1,0))
    def test_physical_bound_has_scale_free_positive_coefficients(self):self.assertEqual(b.physical_bound((1,2,3,4),2,3),3*2*3+2**3*4)
    def test_radial_cutoff_low_orders(self):self.assertEqual(b.radial_cutoff((1,2,5),F(1,4),2),(1,4,22))
    def test_third_radial_cutoff_keeps_all_three_shape_jets(self):self.assertEqual(b.radial_cutoff((1,2,5,7),1,2)[3],F(26,8))
    def test_product_leibniz_cross_terms(self):self.assertEqual(b.product_jets((1,2,5),(7,11,13)),(7,25,92))
    def test_actual_first_three_bump_budgets_are_not_changed(self):self.assertEqual(b.repair_bump_bounds(2),(320,121600,220000000))
    def test_six_explicit_actual_bump_budgets(self):self.assertEqual(len(self.r['actual_repair_bump_log_bounds']),6);self.assertEqual(b.repair_bump_bounds(3)[3],18870643200000)
    def test_logistic_first_and_second_polynomials(self):self.assertEqual(b.logistic_polynomial(1),(0,1,-1));self.assertEqual(b.logistic_polynomial(2),(0,1,-3,2))
    def test_logistic_quotient_constant_term_and_degree(self):
        for k in range(9):
            v=b.logistic_polynomial(k);self.assertEqual(v[0],0);self.assertEqual(v[1],1);self.assertEqual(len(v),k+2)
    def test_generic_Bell_bound_first_order_independent(self):self.assertEqual(b.generic_sigma_bound(1),128)
    def test_integer_partition_multiplicities(self):self.assertEqual(set(b.partitions(4)),{(1,1,1,1),(1,1,2),(1,3),(2,2),(4,)})
    def test_generic_bounds_exist_in_positive_integers_through_eight(self):
        for k in range(9):self.assertIs(type(b.generic_sigma_bound(k)),int);self.assertGreater(b.generic_sigma_bound(k),0)
    def test_old_sharper_sigma_bounds_used_when_available(self):self.assertEqual([b.sigma_bound(k) for k in range(4)],[1,32,1216,220000])
    def test_U_and_F_bump_inverse_radius_costs(self):
        self.assertEqual(b.bump_radial((7,11,17),3,2,'U'),(7,F(11,3),F(28,9)))
        self.assertEqual(b.bump_radial((7,11,17),3,2,'F'),(F(7,2),F(29,12),F(59,24)))
    def test_low_order_repaired_interface_independent_legacy_comparator(self):self.assertTrue(self.r['checks']['low_order_interface_agrees_with_original_first_seminorms'])
    def test_disjoint_repair_does_not_gain_factor_five(self):self.assertEqual(b.extended_jets((1,),(2,),(7,),3,2,1),(23,))
    def test_nonzero_shape_third_jet_is_not_suppressed(self):
        small=b.radial_cutoff((1,2,5,7),1,2);big=b.radial_cutoff((1,2,5,700),1,2)
        self.assertEqual(big[:3],small[:3]);self.assertGreater(big[3],small[3])
    def test_all_six_auxiliary_controls_preserve_low_jet_budget(self):self.assertEqual(len(self.r['auxiliary_controls']),6);self.assertTrue(all(all(v['checks'].values()) for v in self.r['auxiliary_controls']))
    def test_center_jets_are_exact_not_floating_samples(self):
        c=b.control(10**24);a=b.control_constants()['a'];self.assertEqual(c['center_first_difference'],a/10**24);self.assertEqual(c['center_second_difference'],0);self.assertEqual(c['center_third_total'],3712-a*10**24)
    def test_C1_difference_goes_to_zero_and_C2_difference_stays_bounded(self):
        lo,hi=b.control(1),b.control(10**36);self.assertLess(hi['first_difference'],lo['first_difference']);self.assertLess(hi['second_difference'],lo['second_difference']);self.assertGreater(hi['second_difference'],0)
    def test_central_third_jet_can_exceed_arbitrary_fixed_target(self):
        a=b.control_constants()['a'];target=10**50;n=int((F(target)+3713)/a)+1;self.assertLess(b.control(n)['center_third_total'],-target)
    def test_same_affine_budget_accepts_all_controls(self):self.assertTrue(all(v['all_rows_below_common_budget'] for v in self.r['low_order_affine_controls']))
    def test_full_pressure_baseline_present_in_every_absolute_ledger(self):
        for v in self.r['complete_pressure_difference_controls']:self.assertEqual(v['fixed_pressure_baseline'][2],'11/2');self.assertEqual(v['fixed_pressure_baseline'][4],'13/2')
    def test_repair_difference_is_proportional_to_value_difference(self):self.assertTrue(self.r['checks']['repair_difference_is_uniformly_value_small'])
    def test_theta_viscous_primitive_sign_and_measure(self):self.assertTrue(self.r['symbolic']['checks']['theta_viscous_primitive'])
    def test_axial_viscous_primitive_sign_and_measure(self):self.assertTrue(self.r['symbolic']['checks']['axial_viscous_primitive'])
    def test_second_stress_principal_terms_cannot_be_omitted(self):self.assertTrue(self.r['symbolic']['checks']['theta_second_stress_principal']);self.assertTrue(self.r['symbolic']['checks']['axial_second_stress_principal'])
    def test_bad_orders_and_frequencies_rejected(self):
        for v in (-1,True,F(1,2),1.0):
            with self.assertRaises(ValueError):b.log_coefficients(v)
            with self.assertRaises(ValueError):b.physical_coefficients(v)
            with self.assertRaises(ValueError):b.control(v)
        with self.assertRaises(ValueError):b.control(0)
    def test_float_and_bool_bounds_rejected(self):
        for v in (True,.1):
            with self.assertRaises(TypeError):b.radial_cutoff((1,v),1,1)
            with self.assertRaises(TypeError):b.physical_bound((1,2),v,1)
    def test_missing_or_invalid_jet_inputs_rejected(self):
        for L,d,x in (((1,-1),1,1),((0,2),1,1),((1,2),0,1),((1,2),1,0)):
            with self.assertRaises(ValueError):b.radial_cutoff(L,d,x)
        with self.assertRaises(ValueError):b.product_jets((1,2),(1,))
        with self.assertRaises(ValueError):b.physical_bound((1,2),2,2)
        with self.assertRaises(ValueError):b.extended_jets((1,),(1,),(1,2),1,2,1)
        with self.assertRaises(ValueError):b.extended_jets((1,),(1,),(1,),-1,2,1)
    def test_invalid_bump_radius_or_type_rejected(self):
        for p,r,t in ((1,2,'F'),(0,1,'U'),(1,1,'other')):
            with self.assertRaises(ValueError):b.bump_radial((1,2),p,r,t)
    def test_parent_and_frozen_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'f84076cef9550e7bd146493ce9a5a4f7afdc5b3d');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record_reproduced_exactly(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_bindings_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_scope_promotion_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_is_not_false_formal_flag(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_pressure_baseline_or_nondegeneracy_cannot_be_erased(self):
        for key in ('repair_pressure','nondegeneracy'):
            r=deepcopy(self.record);r['protocol'][key]='omitted'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_original_Borel_order_requirements_cannot_be_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(2)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_cli_replay_and_frozen_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
