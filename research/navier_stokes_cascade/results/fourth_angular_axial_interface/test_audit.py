"""Fourth reference-root jets, complete energy, compact gate and scope controls."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import angular_bounds as b

class FourthAngularTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_reciprocal_partial_fractions(self):self.assertTrue(self.r['symbolic']['checks']['reciprocal_partial_fractions_exact'])
    def test_all_four_log_derivatives(self):self.assertTrue(self.r['symbolic']['checks']['logarithm_derivatives_through_four_exact'])
    def test_log_fourth_norm_fits_same_budget(self):self.assertEqual(self.r['budgets']['constants']['log_C4_upper'],'31/6');self.assertTrue(self.r['checks']['log_fourth_norm_fits_same_eight_budget'])
    def test_same_angular_target_and_root_radius(self):
        c=self.r['budgets']['constants'];self.assertEqual(c['angular_target_C4_coefficient'],2**29);self.assertEqual(c['angular_root_C4_coefficient'],2**32)
    def test_same_pointwise_root_uniqueness(self):self.assertTrue(self.r['checks']['same_pointwise_uniqueness_ball'])
    def test_fourth_interpolation_formula(self):self.assertTrue(self.r['symbolic']['checks']['fourth_interpolation_power_all_terms'])
    def test_raw_fourth_interpolation_cap(self):self.assertEqual(self.r['budgets']['constants']['interpolation_fourth_raw'],408);self.assertTrue(self.r['checks']['full_interpolation_C4_below64'])
    def test_interpolation_energy_T64_endpoint(self):self.assertTrue(self.r['checks']['interpolation_energy_at_T64'])
    def test_power_shape_fourth_factor_is_25_not_24(self):self.assertTrue(self.r['checks']['p_squared_factorial_C4_is25'])
    def test_complete_wait_and_tail_energy(self):self.assertTrue(self.r['checks']['unedited_energy_C4_budget']);self.assertTrue(self.r['checks']['full_post_energy_C4_below8192T'])
    def test_both_signed_angular_edits_included(self):self.assertTrue(self.r['checks']['complete_signed_edit_C4_below_one'])
    def test_full_prefix_s0_polynomial_fourth_norm(self):self.assertTrue(self.r['checks']['prefix_s0_polynomial_C4_is93'])
    def test_strength_epsilon_includes_end_squares(self):self.assertTrue(self.r['checks']['complete_coefficient_epsilon_budget'])
    def test_all_four_raw_strength_product_terms(self):self.assertTrue(self.r['symbolic']['checks']['fourth_strength_derivative_all_linear_and_quadratic_terms'])
    def test_linear_coefficient_derivatives_not_removed(self):self.assertTrue(self.r['symbolic']['checks']['dropping_linear_coefficient_derivatives_detected'])
    def test_normalized_twelve_epsilon_caps(self):
        e=F(1,2**20);u=b.strength_jets(e);self.assertEqual(u[0],2);self.assertEqual(u[1],9*e);self.assertTrue(all(v<12*e for v in u[1:]))
    def test_all_four_uniform_caps_are_inductive(self):self.assertTrue(self.r['checks']['uniform_twelve_caps_inductive'])
    def test_zero_epsilon_has_zero_derivative_caps(self):self.assertEqual(b.strength_jets(0),(2,0,0,0,0))
    def test_root_jets_monotone_in_epsilon(self):
        lo=b.strength_jets(F(1,2**22));hi=b.strength_jets(F(1,2**20));self.assertTrue(all(a<=c for a,c in zip(lo,hi)))
    def test_strength_full_C4_under_two(self):self.assertTrue(self.r['checks']['fourth_strength_full_norm_below_two'])
    def test_affine_end_corrections_kept(self):self.assertEqual(self.r['budgets']['constants']['end_coefficient_C4_multiplier'],4)
    def test_later_reference_pulse_C4_is_110(self):self.assertEqual(self.r['budgets']['constants']['reference_axial_pulse_C4_upper'],110)
    def test_all_five_manufactured_fourth_jets(self):
        d=self.r['manufactured'];self.assertEqual(len(d['rows']),5);self.assertTrue(d['manufactured_not_selected_reference_root']);self.assertTrue(d['checks']['all_five_same_root_fourth_jet_controls'])
    def test_manufactured_endpoints_included(self):self.assertEqual([v['eta'] for v in self.r['manufactured']['rows']][::4],['-1','1'])
    def test_fourth_compact_quadratic_root_terms(self):
        c=b.compact_cost(2,3,100,F(1,1000),5,7);self.assertEqual(c['root_fourth_jet'],740);self.assertEqual(c['root_C4'],740+F(1,250))
    def test_compact_known_C3_to_C4_four_jet_conversion(self):
        c=b.compact_cost(1,1,1,1,0,0);self.assertEqual(c['root_C4']-c['root_fourth_jet'],4)
    def test_compact_actual_unit_mass_bump_value(self):self.assertTrue(self.r['checks']['actual_compact_unit_mass_bump_value'])
    def test_compact_primitive_and_transformed_target_different_coefficients(self):
        c=b.compact_cost(2,3,100,0,5,7);self.assertEqual(c['missing_input'],44800005)
    def test_compact_gate_exact_at_boundary(self):
        c=b.compact_cost(2,3,100,F(1,1000),5,7);g=b.compact_gate(2,3,100,F(1,1000),5,7,c['repaired_U_C4']);self.assertEqual(g['lhs'],g['rhs']);self.assertTrue(g['accepts'])
    def test_compact_gate_rejects_less_than_majorant(self):self.assertTrue(self.r['checks']['too_small_compact_certificate_rejected'])
    def test_compact_gate_equivalent_over_rational_inputs(self):
        for A in (1,2):
            for N in (1,7):
                for r in (0,F(1,13)):
                    for K in (0,17,F(10**15)):
                        g=b.compact_gate(A,3,N,r,5,11,K);self.assertEqual(g['accepts'],g['same_positive_envelope_inequality'])
    def test_fourth_target_not_dropped_when_primitive_zero(self):
        c=b.compact_cost(2,3,100,0,0,7);self.assertGreater(c['repaired_U_C4'],3)
    def test_lower_jet_nonlinearity_not_dropped_when_target_zero(self):
        c=b.compact_cost(2,3,100,F(1,1000),0,0);self.assertEqual(c['root_fourth_jet'],600)
    def test_axis_fourth_Cauchy_factorial_cancellation(self):
        a=b.axis_fourth(2,F(1,8),4,F(1,100),F(1,10000));old=b.module('collar_operator_envelope','envelope.py').bounds(2,F(1,8),4,F(1,100),F(1,10000));self.assertEqual(a['factor'],old['radial_factor']*sum(old['delta']**(-j) for j in range(5)))
    def test_global_cover_does_not_stop_at_compact_annulus(self):self.assertEqual(self.r['interface']['pressure_fixture']['Kleading4'],'110')
    def test_full_pressure_weighted_axial_term_retained(self):self.assertTrue(self.r['checks']['both_signed_axial_viscosities_retained'])
    def test_global_first_C3_remains_an_input(self):
        a=b.pressure_interface(F(1,200),7,13,17,19);c=b.pressure_interface(F(1,200),7,13,17,23);self.assertGreater(c['I0'],a['I0']);self.assertGreater(c['I1'],a['I1'])
    def test_abstract_C3_does_not_supply_fourth_jet(self):self.assertTrue(self.r['checks']['abstract_small_C3_does_not_bound_C4'])
    def test_invalid_epsilon(self):
        for e in (-1,F(1,2**19)):
            with self.assertRaises(ValueError):b.strength_jets(e)
    def test_float_bool_inputs_rejected(self):
        for e in (0.,False):
            with self.assertRaises(TypeError):b.strength_jets(e)
    def test_invalid_compact_inputs(self):
        for args in ((0,1,1,0,0,0),(1,0,1,0,0,0),(1,1,0,0,0,0),(1,1,1,-1,0,0),(1,1,1,0,-1,0),(1,1,1,0,0,-1)):
            with self.assertRaises(ValueError):b.compact_cost(*args)
    def test_negative_trial_gate_rejected(self):
        with self.assertRaises(ValueError):b.compact_gate(1,1,1,0,0,0,-1)
    def test_invalid_pressure_global_inputs(self):
        for args in ((F(1,200),7,-1,1,1),(0,7,1,1,1),(F(1,200),0,1,1,1)):
            with self.assertRaises(ValueError):b.pressure_interface(*args)
    def test_parent_protocol_and_frozen_record(self):self.assertEqual(audit.protocol()['parent_commit'],'5ad730f469503d4d817e3abcd6d58390d8281152');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_history_byte_bindings_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_unverified_scope_promotions_rejected(self):
        for key,value in audit.FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_is_not_false(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_heat_energy_pressure_and_target_normalizations_not_erased(self):
        for key in ('heat','energy','same_repair','compact','global_pressure'):
            r=deepcopy(self.record);r['protocol'][key]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_actual_B4_Z4_obligation_not_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_evidence_and_frozen_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
