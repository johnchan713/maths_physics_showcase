"""Exact jet controls and refusal to promote unquantified continuum inputs."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import input_bounds as b

class FourthInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_all_35_principal_rows_exact(self):self.assertTrue(self.r['principal']['checks']['all_35_exact_principal_splits'])
    def test_every_moment_axis_and_combination_checked(self):self.assertEqual({r['mode'] for r in self.r['principal']['rows']},set(range(7)))
    def test_eta_endpoints_center_interior_covered(self):self.assertEqual({r['eta'] for r in self.r['principal']['rows']},{'-1','-1/2','0','2/5','1'})
    def test_fifth_M_axial_coefficient_retains_U(self):
        a=b.principal(2,3,1,Q(2,3),Q(2,5),Q(1,200),[1,0,0,0,0,0]);c=b.principal(2,3,1,0,Q(2,5),Q(1,200),[1,0,0,0,0,0]);self.assertNotEqual(a['p2'],c['p2'])
    def test_factor_five_not_one(self):
        r=b.principal(2,3,1,0,0,Q(1,200),[0,0,0,0,1,0]);self.assertEqual(r['p2'],-10)
    def test_I5_only_affects_p1(self):
        r=b.principal(2,3,1,0,Q(1,2),Q(1,200),[0,1,0,0,0,0]);self.assertNotEqual(r['p1'],0);self.assertEqual(r['p2'],0)
    def test_J5_only_affects_p1(self):
        r=b.principal(2,3,1,0,0,Q(1,200),[0,0,1,0,0,0]);self.assertEqual(r['p1'],-Q(5,3));self.assertEqual(r['p2'],0)
    def test_S5_only_affects_p2(self):
        r=b.principal(2,3,1,0,0,Q(1,200),[0,0,0,1,0,0]);self.assertEqual(r,dict(p1=0,p2=-5))
    def test_Cp5_and_axis5_both_retained(self):self.assertTrue(self.r['principal']['checks']['both_pressure_datum_and_Cp_fifth_rows_retained'])
    def test_endpoints_do_not_erase_all_fifth_inputs(self):self.assertTrue(self.r['principal']['checks']['endpoint_degeneracy_does_not_remove_I5_or_M5'])
    def test_fifth_pressure_rows_cancel_only_as_actual_sum(self):
        self.assertEqual(b.principal(2,3,1,0,0,Q(1,200),[0,0,0,0,3,-3]),dict(p1=0,p2=0))
    def test_lower_denominator_jets_not_frozen(self):self.assertTrue(self.r['principal']['checks']['full_denominator_lower_jets_not_frozen'])
    def test_exact_jet_product_and_inverse(self):
        a=b.Jet([2,3,5,7,11]);self.assertEqual((a*a.inverse()).a,(1,0,0,0,0));self.assertEqual((b.Jet([0,1])**4).a,(0,0,0,0,1))
    def test_normalized_moment_derivative_has_factor_five(self):self.assertEqual(b.derivative([0,1,2,3,4,5]).a,(1,4,9,16,25))
    def test_missing_fifth_data_rejected(self):
        with self.assertRaises(ValueError):b.derivative([0]*5)
    def test_singular_inverse_rejected(self):
        with self.assertRaises(ValueError):b.Jet(0).inverse()
    def test_invalid_compact_inputs_rejected(self):
        for change in (dict(X=0),dict(H0=0),dict(E=[0]),dict(eta=2),dict(h=0)):
            with self.assertRaises(ValueError):b.coordinates(**dict(b.fixture(),**change))
    def test_compatible_symbolic_linear_growth_limit(self):self.assertTrue(self.r['compatible']['checks']['principal_growth_limit_exact'])
    def test_lower_terms_cannot_cancel_linear_growth(self):self.assertTrue(self.r['compatible']['checks']['nonzero_linear_growth_not_cancelled_by_lower_terms'])
    def test_four_exact_compatible_forward_controls(self):self.assertTrue(self.r['compatible']['checks']['four_exact_original_forward_map_controls'])
    def test_all_five_genuine_forward_moments_documented(self):self.assertIn('Cp=log(2)/2',self.r['compatible']['all_five_forward_moments'])
    def test_control_has_uniform_C4_native_bound(self):self.assertEqual(self.r['compatible']['uniform_U_Ceta4_bound'],'(65/24)*eps for n>=1')
    def test_control_is_not_the_actual_profile_or_cone_failure(self):self.assertTrue(self.r['compatible']['manufactured_not_actual_profile_or_cone_family'])
    def test_exact_control_scaled_error_refines(self):self.assertTrue(self.r['compatible']['checks']['normalized_growth_error_decreases'])
    def test_fourth_composition_all_mixed_partitions(self):self.assertTrue(self.r['moving']['checks']['all_fourth_composition_partitions_exact'])
    def test_implicit_inverse_uses_mixed_remainder(self):self.assertTrue(self.r['moving']['checks']['moving_inverse_mixed_remainder_is_nonzero'])
    def test_new_inverse_fourth_value(self):self.assertEqual(self.r['moving']['pulled_back_normalized_fourth'],'135/4096')
    def test_independent_Lagrange_inversion(self):self.assertEqual(self.r['moving']['lagrange_inversion_fourth'],'135/4096')
    def test_prior_explicit_lift_not_confused_with_inverse(self):self.assertEqual(self.r['moving']['previous_explicit_lift_normalized_fourth'],'3/4096')
    def test_internal_phase_freezing_rejected(self):self.assertTrue(self.r['moving']['checks']['freezing_internal_phase_is_rejected'])
    def test_inverse_and_composition_agree(self):self.assertTrue(self.r['moving']['checks']['independent_composition_matches_inverse_recursion'])
    def test_all_inverse_jets_retained(self):self.assertEqual(self.r['moving']['inverse_normalized_jets'],['-sqrt(3)/4','sqrt(3)/16','sqrt(3)/128','-sqrt(3)/64'])
    def test_parent_protocol_and_frozen_evidence(self):self.assertEqual(audit.protocol()['parent_commit'],'187149d4cdd64da1d82d405dc3cfbfd3c92519f5');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_bindings(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_unproved_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_is_not_boolean_false(self):
        r=deepcopy(self.record);r['flags']['unconditional_foundations_blowup_or_formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_pressure_gauge_root_phase_and_gate_not_erased(self):
        for k in ('pressure','principal','root','moving','primitive','gate'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_actual_input_obligation_not_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_false_algebra_gate_rejected(self):
        r=deepcopy(self.record);r['gates']['principal_all_35_exact_principal_splits']=False
        with self.assertRaises(ValueError):audit.validate(r)
    def test_frozen_output_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
