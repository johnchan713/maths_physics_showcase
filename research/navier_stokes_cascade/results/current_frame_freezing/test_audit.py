"""Mixed current-frame and global cone controls, not actual profile/scale acceptance."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import mixed_bounds as b

class MixedFreezingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=audit.generate();cls.r=cls.record['reviews']
        cls.data=b.deficit().constants(F(1,2),2,2,6,4,F(1,2),F(3,4))
        cls.cost=b.current_costs(cls.data,1,20,6,F(1,20000));cls.cert=b.deficit().certificate(cls.data,F(1,100))
    def controls(self,k):self.assertTrue(all(r['checks'][k] for r in self.r['endpoints']['rows']))
    def symbolic(self,k):self.assertTrue(self.r['symbolic']['checks'][k])
    def test_all_192_current_endpoint_controls(self):self.assertEqual(self.r['endpoints']['case_count'],192);self.assertTrue(all(self.r['endpoints']['checks'].values()))
    def test_leading_floor_and_raw_deficit_controls(self):self.controls('leading_endpoint_and_global_floor_controls')
    def test_leading_two_point_jets(self):self.controls('leading_two_point_jets')
    def test_both_current_endpoint_value_errors(self):self.controls('both_current_endpoint_value_errors')
    def test_relative_current_target_error(self):self.controls('pointwise_relative_target_error')
    def test_two_point_current_shear_triangle_bound(self):self.controls('current_two_point_shear_bound')
    def test_current_endpoint_value_spectral_and_mixing_floors(self):self.controls('current_endpoint_floors')
    def test_current_orientation_not_leading_omega(self):self.controls('actual_unit_direction_orientation_spends_omega_over_three')
    def test_unfrozen_same_tilt_margin(self):self.controls('unfrozen_same_tilt_margin')
    def test_relative_V_mixed_bound(self):self.controls('mixed_relative_V_bound')
    def test_relative_W_mixed_bound(self):self.controls('mixed_relative_W_bound')
    def test_combined_error_spends_half_margin_budget(self):self.controls('mixed_freezing_budget_fits_half_margin')
    def test_frozen_half_margin_positive(self):self.controls('frozen_same_tilt_half_margin')
    def test_stress_amplitude_cancels_at_two_magnitudes(self):
        self.controls('target_amplitude_cancels_exactly');self.assertEqual(len({r['target_scale'] for r in self.r['endpoints']['rows']}),2)
    def test_three_original_tilts_and_all_eight_signs(self):
        rows=self.r['endpoints']['rows'];self.assertEqual({r['u'] for r in rows},{'3/4','4/3','5/12'});self.assertEqual(len({tuple(r['signs']) for r in rows}),8)
    def test_both_target_components_retained(self):self.assertTrue(self.r['endpoints']['checks']['pointwise_target_has_both_components'])
    def test_rational_mixing_lower_and_upper(self):
        self.assertEqual(self.cost['mixing_squared_lower'],F(1,21));self.assertEqual(self.cost['mixing_squared_upper'],182)
        self.assertLessEqual(self.cost['endpoint_inputs']['cmin']**2,F(1,21));self.assertGreaterEqual(self.cost['endpoint_inputs']['cmax']**2,182)
    def test_all_current_floor_inputs_are_spent(self):
        p=self.cost['endpoint_inputs'];self.assertEqual((p['f'],p['Fbar'],p['t'],p['Gbar'],p['omega']),(F(1,4),3,1,7,F(1,6)))
    def test_current_slow_mesh_cost_explicit(self):self.assertEqual(self.cost['Cmesh'],71129772)
    def test_perturbation_cost_counts_two_endpoints(self):self.assertEqual(self.cost['Cpert'],3332196);self.assertEqual((self.cost['value_error_cost']['inputs']['JF'],self.cost['value_error_cost']['inputs']['JG']),(2,2))
    def test_each_coordinate_cost_is_linear_before_maximum(self):
        p=self.cost['endpoint_inputs'];args=tuple(p[k] for k in ('f','Fbar','t','Gbar','omega','cmin','cmax','u'))
        total=b.leading().frame_cost(*args,F(3,5)+F(4,7),F(11,3)+F(7,11))
        a=b.leading().frame_cost(*args,F(3,5),F(11,3));c=b.leading().frame_cost(*args,F(4,7),F(7,11))
        for k in ('relative_V_Lipschitz','relative_W_Lipschitz'):self.assertEqual(total[k],a[k]+c[k])
        self.assertLessEqual(total['relative_coordinate_Lipschitz'],a['relative_coordinate_Lipschitz']+c['relative_coordinate_Lipschitz'])
    def test_zero_leading_jets_do_not_erase_current_value_cost(self):
        c=b.current_costs(self.data,0,0,6,F(1,20000));self.assertEqual(c['Cmesh'],0);self.assertGreater(c['Cpert'],0)
    def test_band_normalization_is_eight(self):self.assertEqual(self.cost['profile_normalization'],8);self.symbolic('normalization_eight_below_safe_three_halves_power')
    def test_background_required_deficit_has_sixteen(self):
        E=F(1,10**30);r=b.mixed_gate(self.data,self.cost,E,10000)
        self.assertEqual(r['background_deficit'],16*self.data['Lpsi']*E/self.data['inputs']['omega']**2);self.symbolic('background_deficit_has_sixteen')
    def test_freezing_deficit_includes_mesh_and_value_errors(self):
        E=F(1,10**30);r=b.mixed_gate(self.data,self.cost,E,10000)
        self.assertEqual(r['relative_freezing_error_upper'],self.cost['Cmesh']/F(10000**6)+8*self.cost['Cpert']*E)
        self.assertEqual(r['freezing_deficit'],72*self.data['kappa']*r['relative_freezing_error_upper']/self.data['inputs']['omega']**2);self.symbolic('mixed_freezing_deficit_has_both_terms')
    def test_small_background_error_does_not_erase_mesh(self):
        r=b.mixed_gate(self.data,self.cost,0,2);self.assertEqual(r['perturbation_freezing_error'],0);self.assertGreater(r['mesh_freezing_error'],0)
    def test_zero_mesh_error_does_not_certify_current_freezing(self):
        c=b.current_costs(self.data,0,0,6,F(1,20000));gamma=self.cert['floors']['unfrozen_cone_margin'];E=gamma/(32*c['Cpert'])
        r=b.mixed_gate(self.data,c,E,10000);self.assertLessEqual(r['profile_error'],self.cert['radius'])
        self.assertEqual(r['mesh_freezing_error'],0);self.assertEqual(r['relative_freezing_error_upper'],gamma/4);self.assertGreater(r['relative_freezing_error_upper'],gamma/8)
    def test_enlarged_positive_chart_gate_is_independent(self):
        c=b.current_costs(self.data,0,0,6,F(1,20000));r=b.mixed_gate(self.data,c,0,2)
        self.assertFalse(r['chart_gate']);self.assertTrue(r['background_base_gate']);self.assertEqual(r['required_deficit'],0)
    def test_too_large_error_fails_background_base(self):self.assertFalse(b.mixed_gate(self.data,self.cost,1,10000)['background_base_gate'])
    def test_too_large_freezing_requirement_cannot_fit_raw_supremum(self):self.assertFalse(b.mixed_gate(self.data,self.cost,0,2)['deficit_requirement_within_possible_range'])
    def test_outer_angular_increment_improves_orientation(self):self.symbolic('nonnegative_outer_angular_stress_improves_orientation')
    def test_outer_angular_increment_reduces_same_tilt_ratio(self):self.symbolic('nonnegative_outer_angular_stress_reduces_cone_ratio')
    def test_outer_relative_magnitude_upper_bound_not_inferred(self):self.symbolic('outer_stress_magnitude_ratio_can_be_unbounded');self.assertFalse(audit.FLAGS['global_outer_relative_stress_smallness_or_magnitude_upper_bound_established'])
    def test_rational_square_root_enclosure_safety(self):self.symbolic('rational_mixing_enclosures_are_safe')
    def test_four_value_only_oscillation_controls(self):
        self.assertEqual(len(self.r['obstruction']['rows']),4);self.assertTrue(all(self.r['obstruction']['checks'].values()))
    def test_pure_mesh_slope_diverges_with_small_value_errors(self):self.assertTrue(self.r['obstruction']['checks']['shrinking_value_error_and_distance_but_growing_slope'])
    def test_all_joint_later_band_gates(self):self.assertTrue(all(self.r['join']['checks'].values()))
    def test_later_witness_splits_cone_budget_into_sixteenths(self):
        w=self.r['join']['witness'];self.assertTrue(w['checks']['mesh_spends_at_most_one_sixteenth']);self.assertTrue(w['checks']['value_perturbation_spends_at_most_one_sixteenth']);self.symbolic('split_mesh_and_perturbation_spend_one_eighth')
    def test_full_pressure_and_shifted_viscosity_costs_not_erased(self):self.assertTrue(self.r['join']['checks']['full_pressure_and_viscosity_rows_increase_actual_error_budget'])
    def test_doubled_box_and_current_orientation_represented(self):self.assertTrue(self.r['join']['checks']['current_orientation_is_omega_over_three']);self.assertTrue(self.r['join']['checks']['doubled_box_and_positive_chart_retained'])
    def test_no_current_or_unit_stress_derivative_input(self):self.assertTrue(self.r['join']['checks']['no_current_derivative_or_unit_target_derivative_input'])
    def test_invalid_jet_box_or_chart_inputs(self):
        for a,c,k,d in ((-1,20,6,F(1,20000)),(1,-1,6,F(1,20000)),(1,20,-1,F(1,20000)),(1,20,6,0)):
            with self.assertRaises(ValueError):b.current_costs(self.data,a,c,k,d)
    def test_invalid_error_and_band(self):
        for E,ell in ((-1,10),(0,1),(0,True)):
            with self.assertRaises(ValueError):b.mixed_gate(self.data,self.cost,E,ell)
    def test_exact_inputs_reject_float_and_bool(self):
        for v in (True,.1):
            with self.assertRaises(TypeError):b.current_costs(self.data,v,20,6,F(1,20000))
            with self.assertRaises(TypeError):b.mixed_gate(self.data,self.cost,v,10)
    def test_invalid_raw_deficit(self):
        for v in (0,-1,self.data['kappa']+1):
            with self.assertRaises(ValueError):b.deficit().certificate(self.data,v)
    def test_conditional_witness_rejects_impossible_shape_budgets(self):
        g=b.deficit().gates();p=g.fixture();rows=g.affine_gate(**p)['rows']
        with self.assertRaises(ValueError):b.later_witness(self.data,self.cost,rows,1,4,F(1,100),p['h'])
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'cd07a79f66c18f731a2a634eb35c7ee6ffffbd63');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record_exact_reproduction(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_bindings_preserved(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_every_unsupported_scope_promotion_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_scope_bool_not_replaced_by_zero(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_mixed_error_cost_tampering_rejected(self):
        r=deepcopy(self.record);r['reviews']['join']['costs']['Cpert']='0'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_nonlinear_closeness_obligation_cannot_be_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(2)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_pressure_premise_cannot_be_removed(self):
        r=deepcopy(self.record);r['flags']['full_corrected_width_pressure_viscosity_and_common_velocity_stress_Borel_inputs_retained']=False
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_cli_replay_and_no_frozen_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
