"""Complete source, signed moment boundary controls and strict conditional scope."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import moments as m

class HigherPressureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_original_ordered_density_equals_flux_all_48_cases(self):
        rows=self.r['exact']['rows'];self.assertEqual(len(rows),48);self.assertTrue(all(r['checks']['ordered_pointwise_flux'] for r in rows))
    def test_192_integrated_scalar_jet_identities(self):self.assertEqual(self.r['exact']['scalar_integrated_jet_identities'],192)
    def test_finite_endpoint_identity_retains_axial_boundary_derivatives(self):self.assertTrue(all(r['checks']['finite_endpoint_moment_identity'] for r in self.r['exact']['rows']))
    def test_full_endpoint_boundary_zero(self):self.assertTrue(all(r['checks']['full_endpoint_vanishes'] for r in self.r['exact']['rows']))
    def test_positive_order_axis_zero_is_stronger_than_regularity(self):self.assertTrue(all(r['checks']['positive_order_zero_axis_trace'] for r in self.r['exact']['rows']))
    def test_regular_nonzero_axis_flux_control(self):
        z=m.identities(1,1,F(0),F(1,200),1,1);self.assertFalse(m.zero(z['axis']))
    def test_actual_positive_order_trace_removes_axis_flux_only(self):
        z=m.identities(1,0,F(0),F(1,200),1,1);self.assertTrue(m.zero(z['axis']));self.assertTrue(z['shifted_omission_detected'])
    def test_pressure_order_two_shift_not_erased(self):self.assertTrue(self.r['exact']['checks']['order_two_pressure_shift_is_nonzero'])
    def test_truncated_future_support_control(self):self.assertTrue(self.r['exact']['checks']['future_support_cannot_be_truncated'])
    def test_both_angular_endpoints_no_division_by_d(self):self.assertTrue(self.r['exact']['checks']['angular_endpoints_included'])
    def test_general_order_seven_not_only_stored_fixtures(self):
        z=m.identities(7,0,F(2,5),F(1,100),F(3,2),F(1,2));self.assertTrue(z['pointwise']);self.assertTrue(z['exact'])
    def test_original_hierarchy_27_independent_float_samples(self):self.assertEqual(self.r['original_hierarchy']['samples'],27);self.assertTrue(self.r['original_hierarchy']['checks']['original_hierarchy_density_matches'])
    def test_original_hierarchy_shift_omission_control(self):self.assertTrue(self.r['original_hierarchy']['checks']['original_axial_viscosity_omission_detected'])
    def test_k0_bound_exactly_reduces_to_old_formula(self):self.assertTrue(self.r['norms']['checks']['k0_exactly_recovers_frozen_signed_bound'])
    def test_ordered_B_and_C_cross_terms(self):
        z=m.source_bounds(1,F(1,200),1,7,(2,3),(5,11),17);self.assertEqual((z['B2'],z['C2']),(37,110))
    def test_weighted_radial_cancellation_does_not_remove_axial_term(self):self.assertTrue(self.r['norms']['checks']['weighted_axial_viscosity_not_radial_boundary_cancellation'])
    def test_unweighted_axis_and_axial_costs(self):self.assertTrue(self.r['norms']['checks']['unweighted_axis_and_axial_terms_both_charged'])
    def test_weighted_nu_difference_is_exact_axial_cap(self):self.assertTrue(self.r['norms']['checks']['weighted_difference_exactly_axial'])
    def test_zero_previous_input_removes_only_axial_cap(self):
        z=m.source_bounds(1,F(1,200),1,7,(2,3),(5,11),0);self.assertEqual(z['axial'],(0,0));self.assertGreater(z['I0'],0)
    def test_missing_previous_angular_input_cannot_default_to_zero(self):
        with self.assertRaises(ValueError):m.source_bounds(1,F(1,200),1,1,(1,2),(1,2))
    def test_large_order_absolute_lambda_difference(self):
        z=m.source_bounds(100,F(1,100),1,1,(F(1),)*101,(F(1),)*101,1);self.assertEqual(z['B2'],101);self.assertGreater(z['I1'],0)
    def test_global_v3_requires_U4_not_U3_in_average_bound(self):self.assertTrue(self.r['norms']['checks']['global_v3_from_U4_retains_extra_angular_jet'])
    def test_reciprocal_third_bound_over_full_h_interval(self):self.assertTrue(self.r['norms']['checks']['reciprocal_third_angular_bound_less_than_two'])
    def test_all_four_manufactured_fixture_caps(self):self.assertEqual(len(self.r['fixture_caps']['rows']),4);self.assertTrue(self.r['fixture_caps']['checks']['all_fixture_caps'])
    def test_second_repair_keeps_three_lower_products(self):
        z=m.second_repair_bounds(2,3,5,7,11,13,17,19,23)
        self.assertEqual(z,(10,20,F(1313,2),2349,1999))
    def test_second_repair_F1_square_not_current_quadratic(self):self.assertTrue(self.r['norms']['checks']['lower_F1_square_retained_in_Cp'])
    def test_second_repair_U1_square_retained(self):self.assertTrue(self.r['norms']['checks']['lower_U1_square_retained_in_S'])
    def test_original_inverse_800_normalization_and_slope_cost(self):self.assertTrue(self.r['norms']['checks']['unchanged_linear_inverse_800_with_amplitude_and_slope_cost'])
    def test_wrong_array_lengths_or_negative_bounds_rejected(self):
        for u,v,p in (((1,),(1,2),0),((1,2),(1,),0),((1,-1),(1,2),0),((1,2),(1,2),-1)):
            with self.assertRaises(ValueError):m.source_bounds(1,F(1,200),1,1,u,v,p)
    def test_invalid_parameters_and_orders_rejected(self):
        for k,h,nu,x in ((-1,F(1,200),1,1),(True,F(1,200),1,1),(1,0,1,1),(1,F(1,50),1,1),(1,F(1,200),-1,1),(1,F(1,200),1,0)):
            with self.assertRaises(ValueError):m.source_bounds(k,h,nu,x,(1,2),(1,2),1)
    def test_no_float_or_bool_norm_inputs(self):
        for bad in (.5,True):
            with self.assertRaises(TypeError):m.source_bounds(1,F(1,200),1,1,(1,bad),(1,2),1)
    def test_regular_average_invalid_inputs(self):
        for r,h,U in ((1,F(1,200),1),(4,F(1,200),1),(2.,F(1,200),1),(2,0,1),(2,F(1,200),-1)):
            with self.assertRaises(ValueError):m.regular_v_bound(0,h,r,U)
    def test_repair_negative_norm_or_zero_support_rejected(self):
        with self.assertRaises(ValueError):m.second_repair_bounds(0,1,1,1,1,1,1,1,1)
        with self.assertRaises(ValueError):m.second_repair_bounds(1,1,1,1,1,1,1,-1,1)
    def test_parent_protocol_and_frozen_record(self):self.assertEqual(audit.protocol()['parent_commit'],'39dca194ea7dbe3f8e505ecaaa42729ca5e64754');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_source_bindings_preserved(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_unsupported_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_cannot_replace_false(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_pressure_support_and_angular_input_protocol_not_erased(self):
        for k in ('support','norms','hierarchy'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_H23_and_global_input_obligations_retained(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_no_frozen_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
