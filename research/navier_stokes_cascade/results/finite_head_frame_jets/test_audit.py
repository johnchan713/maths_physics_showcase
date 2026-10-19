"""Original diagonal quantifiers, physical frame jets and strict conditional scope."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import frame_jets as b

class FiniteHeadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_order_one_recovers_original_half_tail(self):self.assertEqual(b.borel_majorant(1,(7,),1)['q2h_normalized'],F(15,2))
    def test_order_two_recovers_original_half_tail(self):self.assertEqual(b.borel_majorant(2,(7,),1)['q2h_normalized'],F(15,2))
    def test_order_three_has_two_head_coefficients(self):
        r=b.borel_majorant(3,(2,3),F(1,2));self.assertEqual(r['first_tail_order'],3);self.assertEqual(r['head_terms'],(F(1,2),F(3,16)));self.assertEqual(r['tail'],F(1,48))
    def test_order_three_exact_normalized_formula(self):self.assertEqual(b.borel_majorant(3,(2,3),F(1,2))['q2h_normalized'],F(17,6))
    def test_order_three_whole_domain_constant_is_quarter_tail(self):self.assertEqual(b.borel_majorant(3,(2,3),1)['whole_domain_constant'],F(21,4))
    def test_order_four_has_three_finite_head_coefficients(self):self.assertEqual(b.borel_majorant(4,(1,2,3),1)['q2h_normalized'],F(49,8))
    def test_all_fixed_orders_have_correct_finite_head_and_tail(self):
        for m in range(1,9):
            J=max(2,m);r=b.borel_majorant(m,tuple(range(1,J)),1)
            self.assertEqual(r['q2h_normalized'],r['whole_domain_constant']);self.assertEqual(len(r['head_terms']),J-1)
    def test_majorant_monotone_on_whole_small_q_domain(self):
        vals=[b.borel_majorant(3,(2,3),z)['q2h_normalized'] for z in (F(1,100),F(1,10),F(1,2),F(1))]
        self.assertEqual(vals,sorted(vals));self.assertEqual(len(set(vals)),4)
    def test_geometric_tail_independent_partial_sums(self):
        z=F(1,3);r=b.borel_majorant(3,(0,0),z);partial=sum((z/2)**n for n in range(3,20));self.assertLess(partial,r['tail']);self.assertEqual(r['tail']-partial,(z/2)**20/(1-z/2))
    def test_second_coefficient_cost_not_silently_absorbed_into_half(self):self.assertGreater(b.borel_majorant(3,(0,100),1)['q2h_normalized'],F(1,2))
    def test_all_four_profile_controls_pass_lower_diagonal_budget(self):
        rows=self.r['second_coefficient_profile_controls'];self.assertEqual(len(rows),4);self.assertTrue(all(all(v['checks'].values()) for v in rows))
    def test_naive_third_jet_bound_fails_for_large_second_coefficient(self):self.assertTrue(self.r['second_coefficient_profile_controls'][-1]['naive_bound_fails'])
    def test_finite_head_restores_valid_bound(self):
        r=b.second_coefficient_control(10**30);self.assertGreater(r['third_perturbation_at_center'],r['naive_first_coefficient_plus_half_bound']);self.assertLessEqual(r['third_perturbation_at_center'],r['corrected_bound'])
    def test_low_order_profile_budget_independent_of_oscillation_frequency(self):self.assertEqual(b.second_coefficient_control(1)['low_order_stage_norm'],b.second_coefficient_control(10**30)['low_order_stage_norm'])
    def test_radius_derivative_budgets(self):self.assertEqual(b.radius_jets(1,4),((4,1,1,3),(1,1,3,15)))
    def test_first_E_coefficient_uses_full_radius_product(self):self.assertEqual(b.first_coefficient_norm((2,3,5,7),1,4)['E'],(8,14,28,58))
    def test_native_inverse_radius_third_jet(self):self.assertEqual(b.profile_error_jets(1,1,4)['F'],(1,2,6,28))
    def test_lower_shear_terms_are_not_discarded(self):self.assertEqual(b.profile_error_jets(1,1,4)['g_theta'],(32,100,472))
    def test_axial_shear_radius_product_terms(self):self.assertEqual(b.profile_error_jets(1,1,4)['g_z'],(4,5,7))
    def test_Euclidean_shear_error_uses_component_sum(self):self.assertEqual(b.profile_error_jets(1,1,4)['g'],(36,105,479))
    def test_independent_symbolic_shear_and_radius_identities(self):
        r=self.r['symbolic']['checks']
        for k in list(r):
            if 'radius' in k or 'shear' in k:self.assertTrue(r[k])
    def test_normal_second_projected_identity(self):self.assertTrue(self.r['symbolic']['checks']['normal_second_projected_identity'])
    def test_second_quotient_and_root_identities(self):self.assertTrue(self.r['symbolic']['checks']['squared_c_second_quotient']);self.assertTrue(self.r['symbolic']['checks']['square_root_second'])
    def test_frame_normal_derivative_caps(self):
        r=b.spectral_bounds((1,2,4,6,4),(3,5),(7,11));self.assertEqual(r['N'],(F(7,4),F(191,16)));self.assertEqual(r['K'],r['N'])
    def test_rational_c_magnitude_enclosures(self):
        r=b.spectral_bounds((1,2,4,6,4),(3,5),(7,11));self.assertEqual((r['cmin'],r['cmax']),(F(1,6),F(15,2)));self.assertLessEqual(r['cmin']**2,F(4,24));self.assertGreaterEqual(r['cmax']**2,F(60,8))
    def test_discriminant_and_P_derivative_terms(self):
        r=b.spectral_bounds((1,2,4,6,4),(3,5),(7,11));self.assertEqual((r['D1'],r['D2'],r['P1'],r['P2']),(148,418,64,188))
    def test_c_second_derivative_keeps_inverse_cube_cost(self):
        r=b.spectral_bounds((1,2,4,6,4),(3,5),(7,11));q1,q2=r['squared_c_derivatives'];self.assertGreater(r['c'][1],q2/(2*r['cmin']));self.assertEqual(r['c'][1]-q2/(2*r['cmin']),q1*q1/(4*r['cmin']**3))
    def test_growth_rate_product_includes_all_three_factors(self):
        r=b.spectral_bounds((1,2,4,6,4),(3,5),(7,11));self.assertEqual(r['growth'],(F(2079,2),F(10844589,8)))
    def test_all_five_manufactured_spectral_profiles(self):self.assertEqual(len(self.r['spectral_diagnostics']['rows']),5);self.assertTrue(all(all(v['checks'].values()) for v in self.r['spectral_diagnostics']['rows']))
    def test_physical_radial_scale_factors(self):self.assertEqual(b.radial_physical((2,3),4),(8,50));self.assertTrue(self.r['symbolic']['checks']['physical_second_radial_scale'])
    def test_physical_factors_match_original_general_chain(self):self.assertTrue(self.r['checks']['physical_factors_match_independent_general_chain'])
    def test_invalid_order_head_and_q_ranges_rejected(self):
        for m,head,z in ((0,(),1),(True,(),1),(3,(1,),1),(3,(1,2),0),(3,(1,2),2),(3,(1,-1),1)):
            with self.assertRaises(ValueError):b.borel_majorant(m,head,z)
    def test_invalid_profile_and_spectral_inputs_rejected(self):
        for lo,hi in ((0,1),(2,1)):
            with self.assertRaises(ValueError):b.radius_jets(lo,hi)
        with self.assertRaises(ValueError):b.profile_error_jets(-1,1,4)
        with self.assertRaises(ValueError):b.spectral_bounds((1,2,4,6,1000),(3,5),(7,11))
        with self.assertRaises(ValueError):b.spectral_bounds((1,2,0,6,4),(3,5),(7,11))
        with self.assertRaises(ValueError):b.spectral_bounds((1,2,4,6,4),(3,),(7,11))
    def test_exact_arithmetic_rejects_float_and_bool(self):
        for v in (.1,True):
            with self.assertRaises(TypeError):b.profile_error_jets(v,1,4)
            with self.assertRaises(TypeError):b.borel_majorant(3,(1,2),v)
    def test_control_and_physical_bad_inputs_rejected(self):
        for n in (0,-1,True):
            with self.assertRaises(ValueError):b.second_coefficient_control(n)
        with self.assertRaises(ValueError):b.second_coefficient_control(1,1)
        with self.assertRaises(ValueError):b.radial_physical((1,2),0)
    def test_parent_protocol_and_frozen_record(self):self.assertEqual(audit.protocol()['parent_commit'],'79e54e2be7c0b250516e1e0cc2eb9462ffb3bf98');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()));audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_sources_preserved(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_unsupported_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_not_false_formal_flag(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_finite_head_or_pressure_cannot_be_erased(self):
        for k in ('finite_head','pressure'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_H23_obligation_cannot_be_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_frozen_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
