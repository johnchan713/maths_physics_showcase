"""Exact scalar, chart and regional controls; no actual coefficient suprema computed."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import weighted_gate as w

class RegionalWeightedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.s=cls.record['reviews']['shared_Borel_tail_and_finite_seminorms']
    def test_exact_geometric_identity(self):self.assertTrue(all(r['checks']['exact_geometric_identity'] for r in self.s['geometric_rows']))
    def test_finite_partial_sum_plus_exact_remainder(self):self.assertTrue(all(r['checks']['finite_plus_remainder'] for r in self.s['geometric_rows']))
    def test_half_constant_on_whole_unit_q_interval(self):self.assertTrue(all(r['checks']['whole_small_q_half_power'] for r in self.s['geometric_rows']))
    def test_half_constant_attained_by_majorant_at_z_one(self):self.assertEqual(w.tail(1),F(1,2))
    def test_small_z_tail_retains_original_second_power(self):self.assertGreaterEqual(w.tail(F(1,10**10)),F(1,4*10**20))
    def test_higher_start_tail_is_stronger(self):self.assertLess(w.tail(F(1,2),5),w.tail(F(1,2),2))
    def test_shared_schedule_retains_both_budgets(self):self.assertTrue(all(r['checks']['both_budgets_included'] for r in self.s['shared_dyadic_schedules']))
    def test_entire_support_has_strict_logarithmic_margin(self):self.assertTrue(all(r['checks']['strict_full_support_margin'] for r in self.s['shared_dyadic_schedules']))
    def test_schedule_geometric_growth(self):self.assertTrue(all(r['checks']['geometric_growth'] for r in self.s['shared_dyadic_schedules']))
    def test_large_stress_seminorm_changes_selected_subbudget(self):
        low=w.schedule(2,F(1,200),1,1);high=w.schedule(2,F(1,200),1,2**100)
        self.assertGreater(high['dyadic_exponent'],low['dyadic_exponent'])
    def test_old_requirements_must_not_be_reduced(self):self.assertGreaterEqual(w.schedule(2,F(1,200),1,1,10000)['dyadic_exponent'],10001)
    def test_independent_stress_budget_control_has_actual_cutoff_one(self):
        r=self.s['independent_stress_budget_control'];self.assertTrue(r['cutoff_equals_one']);self.assertTrue(r['padded_velocity_support_budget_passes'])
    def test_independent_stress_requirement_is_not_optional(self):
        r=self.s['independent_stress_budget_control'];self.assertTrue(r['stress_support_budget_fails']);self.assertTrue(r['violates_velocity_only_inference'])
    def test_componentwise_vector_norm_conversion_is_charged(self):self.assertTrue(self.s['checks']['componentwise_vector_norm_pad_retained'])
    def test_uniform_five_chart_factor_has_exact_integer_certificate(self):self.assertTrue(self.s['checks']['uniform_chart_five_factor'])
    def test_first_coefficient_constants_include_half_tail(self):
        r=w.value_constant(2,3,4,5,1,2);self.assertEqual(r['E'],F(5,2));self.assertEqual(r['W'],F(11,2))
    def test_stress_and_velocity_constants_remain_independent(self):
        a=w.value_constant(0,0,0,1,1,2);b=w.value_constant(0,0,0,10**10,1,2)
        self.assertEqual(a['F'],b['F']);self.assertGreater(b['Cval'],a['Cval'])
    def test_inverse_inner_radius_cost_is_retained(self):self.assertGreater(w.value_constant(2,3,4,5,F(1,100),2)['Cval'],w.value_constant(2,3,4,5,1,2)['Cval'])
    def test_outer_compact_radius_cost_is_retained(self):self.assertGreater(w.value_constant(2,3,4,5,1,100)['Cval'],w.value_constant(2,3,4,5,1,2)['Cval'])
    def chart(self,key):self.assertTrue(all(r['checks'][key] for r in self.record['reviews']['exact_band_radial_pullback']['rows']))
    def test_exact_F_profile_pullback(self):self.chart('F_profile_identity')
    def test_exact_theta_radial_chain_rule(self):self.chart('theta_radial_chain_rule')
    def test_exact_axial_radial_chain_rule(self):self.chart('axial_radial_chain_rule')
    def test_nonunit_chart_factor_not_dropped(self):self.chart('nonunit_chart_factor_retained')
    def test_F_value_constant(self):self.chart('F_constant')
    def test_both_shear_component_constants(self):self.chart('theta_constant');self.chart('axial_constant')
    def test_Euclidean_shear_bound_not_componentwise_misidentified(self):self.chart('Euclidean_constant_sufficient')
    def test_omitted_inverse_radius_derivative_is_detected(self):self.chart('missing_inverse_radius_derivative_detected')
    def test_all_twenty_seven_chart_controls(self):self.assertEqual(len(self.record['reviews']['exact_band_radial_pullback']['rows']),27)
    def test_all_108_outer_controls_keep_same_tilt(self):
        r=self.record['reviews']['outer_monotonic_cover'];self.assertEqual(len(r['rows']),108);self.assertTrue(all(r['checks'].values()))
    def test_large_outer_error_can_improve_cone(self):
        r=self.record['reviews']['outer_monotonic_cover']['large_edge_control'];self.assertTrue(r['sufficient_global_small_error_test_fails']);self.assertTrue(r['nevertheless_same_tilt_cone_improves'])
    def test_nonnegative_cutoff_is_essential(self):self.assertTrue(self.record['reviews']['outer_monotonic_cover']['checks']['nonnegative_cutoff_requirement_is_essential'])
    def test_composed_value_constant_meets_parent_radius(self):self.assertTrue(self.record['reviews']['composed_local_input_to_parent_gap']['checks']['combined_value_gate'])
    def test_composed_exact_norms_match_parent_inputs(self):self.assertTrue(self.record['reviews']['composed_local_input_to_parent_gap']['checks']['exact_input_norms_fit_gate'])
    def test_composed_shear_floors(self):
        r=self.record['reviews']['composed_local_input_to_parent_gap']['checks'];self.assertTrue(r['actual_spectral_floor']);self.assertTrue(r['actual_mixing_floor'])
    def test_composed_oriented_same_tilt_cone(self):
        r=self.record['reviews']['composed_local_input_to_parent_gap']['checks'];self.assertTrue(r['positive_oriented_unit_floor']);self.assertTrue(r['same_tilt_margin'])
    def test_composed_band_keeps_four_Q_cost(self):
        r=self.record['reviews']['composed_local_input_to_parent_gap'];self.assertEqual(r['full_chart_band_exponent'],r['q_point_exponent']+2)
    def test_previous_diagnostic_grid_point_not_sufficiently_certified(self):self.assertTrue(self.record['reviews']['composed_local_input_to_parent_gap']['checks']['previous_grid_point_not_certified'])
    def test_invalid_tail_range_or_start_rejected(self):
        for args in ((0,2),(2,2),(F(1,2),1),(F(1,2),True)):
            with self.assertRaises(ValueError):w.tail(*args)
    def test_invalid_radius_or_negative_seminorm_rejected(self):
        for args in ((-1,0,0,0,1,2),(0,0,0,0,0,2),(0,0,0,0,3,2)):
            with self.assertRaises(ValueError):w.value_constant(*args)
    def test_zero_or_out_of_retained_h_range_rejected(self):
        for h in (0,F(1,2)):
            with self.assertRaises(ValueError):w.schedule(2,h,1,1)
    def test_float_and_bool_not_exact_inputs(self):
        for x in (1.,True):
            with self.assertRaises(TypeError):w.value_constant(x,1,1,1,1,2)
    def test_protocol_and_parent(self):
        self.assertEqual(audit.protocol()['parent_commit'],'0c7ff52961d301880d0078d462ea8fe402328855');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_exact_record_reproduced(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_sources_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_unsupported_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_exact_constant_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['shared_Borel_tail_and_finite_seminorms']['finite_seminorm_constants']['Cval']='1'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_actual_seminorm_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_bool_integer_scope_tamper_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_cli_replay_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
