"""Exact leading-frame freezing controls; not sampled actual profiles or continuum proof."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import freezing_bounds as b

class FreezingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=audit.generate();cls.r=cls.record['reviews'];cls.chart=b.chart_caps(1,4)
        cls.cost=b.frame_cost(F(1,2),2,2,6,F(1,2),F(1,2),1,F(3,4),1,20)
    def control(self,k):self.assertTrue(all(r['checks'][k] for r in self.r['exact']['rows']))
    def symbolic(self,k):self.assertTrue(self.r['symbolic']['checks'][k])
    def test_all_216_controls(self):self.assertEqual(self.r['exact']['case_count'],216);self.assertTrue(all(self.r['exact']['checks'].values()))
    def test_three_original_tilt_choices_and_exact_roots(self):self.assertEqual({r['u'] for r in self.r['exact']['rows']},{'3/4','4/3','5/12'});self.control('exact_tilt_square_root')
    def test_both_target_components_nonzero(self):self.assertTrue(self.r['exact']['checks']['both_pressure_sensitive_target_components_retained'])
    def test_unit_target_and_orientation(self):self.control('unit_target_and_orientation')
    def test_endpoint_value_floors(self):self.control('endpoint_value_and_mixing_floors')
    def test_raw_shear_difference_budget(self):self.control('input_shear_difference_budget')
    def test_exact_mixing_identity(self):self.control('exact_mixing_parameter_identity')
    def test_squared_mixing_difference(self):self.control('squared_mixing_difference_bound')
    def test_normalized_frame_difference(self):self.control('normalized_frame_difference_bound')
    def test_relative_V_coordinate_bound(self):self.control('relative_V_freezing_bound')
    def test_relative_W_coordinate_bound(self):self.control('relative_W_freezing_bound')
    def test_amplitude_cancels_exactly_at_three_magnitudes(self):self.control('target_magnitude_cancels_exactly');self.assertEqual(len({r['target_magnitude'] for r in self.r['exact']['rows']}),3)
    def test_all_eight_sign_choices(self):self.assertEqual(len({tuple(r['signs']) for r in self.r['exact']['rows']}),8)
    def test_three_slow_distances(self):self.assertEqual(len({r['distance'] for r in self.r['exact']['rows']}),3)
    def test_discriminant_product_and_quotient_costs(self):
        self.assertEqual(self.cost['Dbar'],60);self.assertEqual(self.cost['D_difference_Lipschitz'],332)
        self.assertEqual(self.cost['product_difference_Lipschitz'],46);self.assertEqual(self.cost['squared_mixing_Lipschitz'],1546)
    def test_unit_normal_cost(self):self.assertEqual(self.cost['unit_frame_Lipschitz'],20)
    def test_relative_coordinate_costs(self):self.assertEqual(self.cost['relative_V_Lipschitz'],3172);self.assertEqual(self.cost['relative_W_Lipschitz'],F(280,3))
    def test_no_direction_derivative_or_stress_magnitude_input(self):
        self.assertEqual(set(self.cost['inputs']),{'f','Fbar','t','Gbar','omega','cmin','cmax','u','JF','JG'})
    def test_zero_raw_jets_give_zero_frame_cost(self):
        z=b.frame_cost(F(1,2),2,2,6,F(1,2),F(1,2),1,F(3,4),0,0);self.assertEqual(z['relative_coordinate_Lipschitz'],0)
    def test_small_orientation_increases_both_coordinate_costs(self):
        z=b.frame_cost(F(1,2),2,2,6,F(1,4),F(1,2),1,F(3,4),1,20)
        self.assertGreater(z['relative_V_Lipschitz'],self.cost['relative_V_Lipschitz']);self.assertGreater(z['relative_W_Lipschitz'],self.cost['relative_W_Lipschitz'])
    def test_eta_Z_pole_cancellation(self):self.symbolic('exact_eta_Z_cancellation_at_poles')
    def test_angular_shear_native_derivative(self):self.symbolic('native_angular_shear_X_derivative')
    def test_axial_inverse_radius_native_derivative(self):self.symbolic('native_axial_shear_X_inverse_radius_retained')
    def test_exact_relative_V_algebra(self):self.symbolic('covariance_relative_V_exact_difference')
    def test_amplitude_and_fixed_root_cancel(self):self.symbolic('amplitude_and_fixed_tilt_factor_cancel_in_relative_V')
    def test_flat_amplitude_wrong_freezing_countercontrol(self):self.symbolic('flat_amplitude_freezing_exponent_exact');self.symbolic('flat_amplitude_error_does_not_tend_to_zero')
    def test_flat_countercontrol_exponent_stays_above_one(self):
        for y in (F(1,2),F(1,10),F(1,10**20)):
            self.assertGreater((2+y*y)/(1+y*y)**2,1)
    def test_enlarged_profile_rectangle(self):self.assertEqual((self.chart['Xlo'],self.chart['Xhi'],self.chart['Rlo'],self.chart['Rhi']),(F(1,2),32,1,8))
    def test_coordinate_caps(self):self.assertEqual(self.chart['X_coordinate_Lipschitz'],4640);self.assertEqual(self.chart['eta_coordinate_Lipschitz'],14)
    def test_bootstrap_distance_controls_both_exits(self):
        d=self.chart['admissible_slow_distance'];self.assertLessEqual(18*d,F(1,16));self.assertLessEqual(4640*d,F(1,4))
    def test_all_geometry_caps(self):self.assertTrue(all(self.r['geometry']['checks'].values()))
    def test_retained_original_doubled_and_larger_boxes(self):self.assertEqual([b.retained_box_cost(k) for k in (1,2,3)],[4,6,8])
    def test_full_profile_jet_conversion_includes_inverse_radius(self):
        jets=b.profile_to_slow_jets(self.chart,dict(F_C3=1,F_X_C2=2,F_XX_C1=3),dict(U_X_C2=5,U_XX_C1=7))
        self.assertEqual(jets['profile_g_X'],4+192+5+56);self.assertEqual(jets['profile_g_eta'],128+40)
        self.assertEqual(jets['slow_F_Lipschitz'],2*4640+14);self.assertEqual(jets['slow_g_Lipschitz'],257*4640+168*14)
    def test_enlarged_radius_changes_raw_shear_conversion(self):
        j=dict(F_C3=1,F_X_C2=2,F_XX_C1=3);u=dict(U_X_C2=5,U_XX_C1=7)
        self.assertGreater(b.profile_to_slow_jets(b.chart_caps(1,8),j,u)['slow_g_Lipschitz'],b.profile_to_slow_jets(self.chart,j,u)['slow_g_Lipschitz'])
    def test_Cfreeze_and_independent_geometry_budget(self):
        B=b.freezing_budget(self.chart,self.cost,6);self.assertEqual(B['Cfreeze'],6*3172);self.assertEqual(B['chart_distance_threshold'],6/self.chart['admissible_slow_distance'])
    def test_later_freezing_witness_checks_both_gates(self):
        r=b.later_freezing_band(b.freezing_budget(self.chart,self.cost,6),F(1,100));self.assertTrue(all(r['checks'].values()));self.assertEqual(r['frozen_margin'],F(1,200))
    def test_geometry_gate_cannot_be_erased_by_zero_frame_cost(self):
        B=dict(Cfreeze=F(0),chart_distance_threshold=F(10**60));r=b.later_freezing_band(B,1)
        self.assertGreater(r['geometry_ell'],r['freezing_ell']);self.assertTrue(all(r['checks'].values()))
    def test_freezing_gate_can_dominate_geometry_gate(self):
        B=dict(Cfreeze=F(10**60),chart_distance_threshold=F(1));r=b.later_freezing_band(B,F(1,100))
        self.assertGreater(r['freezing_ell'],r['geometry_ell']);self.assertTrue(all(r['checks'].values()))
    def test_complete_global_jet_and_deficit_join(self):self.assertTrue(all(self.r['join']['checks'].values()));self.assertEqual(self.r['join']['budget']['box_diameter_constant'],'6')
    def test_invalid_chart_radii(self):
        for a,c in ((0,1),(-1,4),(2,1)):
            with self.assertRaises(ValueError):b.chart_caps(a,c)
    def test_invalid_frame_value_inputs(self):
        args=[F(1,2),2,2,6,F(1,2),F(1,2),1,F(3,4),1,20]
        for index,value in ((0,0),(1,0),(2,0),(3,1),(4,0),(4,2),(5,0),(6,F(1,4)),(7,0),(8,-1),(9,-1)):
            changed=args[:];changed[index]=value
            with self.assertRaises(ValueError):b.frame_cost(*changed)
    def test_invalid_freezing_distance_enlargement_and_margin(self):
        with self.assertRaises(ValueError):b.freezing_budget(self.chart,self.cost,-1)
        with self.assertRaises(ValueError):b.retained_box_cost(F(1,2))
        for gamma in (0,-1,2):
            with self.assertRaises(ValueError):b.later_freezing_band(dict(Cfreeze=1,chart_distance_threshold=1),gamma)
    def test_invalid_negative_raw_jet_input(self):
        with self.assertRaises(ValueError):b.profile_to_slow_jets(self.chart,dict(F_C3=1,F_X_C2=-1,F_XX_C1=2),dict(U_X_C2=1,U_XX_C1=1))
    def test_floats_and_bools_rejected(self):
        for v in (.1,True):
            with self.assertRaises(TypeError):b.chart_caps(v,4)
            with self.assertRaises(TypeError):b.retained_box_cost(v)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'b0f179b085aa0a64f281fd9ada47951c28e0e196');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record_exact_replay(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_all_historical_bindings_preserved(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_false_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_false_and_zero_types_not_interchangeable(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_nonlinear_frame_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(1)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_Cfreeze_tampering_rejected(self):
        r=deepcopy(self.record);r['reviews']['join']['budget']['Cfreeze']='1'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_erased_pressure_orientation_premise_rejected(self):
        r=deepcopy(self.record);r['flags']['pointwise_full_pressure_sensitive_orientation_and_cone_remain_required']=False
        with self.assertRaises(ValueError):audit.validate(r)
    def test_cli_replay_and_no_frozen_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
