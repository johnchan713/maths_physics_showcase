"""Finite input estimates, exact identities and explicit unsupported-scope guards."""
from copy import deepcopy
from fractions import Fraction as F
import importlib.util,json,tempfile,unittest
from pathlib import Path
import audit
import seminorms as w

class FirstSeminormTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_Cauchy_full_first_angular_norm(self):self.assertEqual(w.collar_jets(2,3,F(1,4),1)['value'],15)
    def test_diagnostic_Picard_upper_bound_and_remaining_strip(self):self.assertTrue(all(self.r['illustrative_Picard_enclosure']['checks'].values()))
    def test_second_angular_derivative_not_first_derivative_budget(self):self.assertEqual(w.collar_jets(2,3,F(1,4),1)['angular_derivative_norm'],108)
    def test_system_rhs_includes_both_matrices_and_source(self):self.assertEqual(w.collar_jets(2,3,F(1,4),1)['system_rhs_norm'],1240)
    def test_F_regular_second_radial_identity(self):self.assertTrue(self.r['exact_radial_and_E_to_F_identities']['F_regular_system'])
    def test_U_regular_second_radial_identity(self):self.assertTrue(self.r['exact_radial_and_E_to_F_identities']['U_regular_system'])
    def test_radial_singular_costs_are_kept_on_positive_radius(self):self.assertGreater(w.collar_jets(2,3,F(1,4),F(1,10))['radial2'],w.collar_jets(2,3,F(1,4),1)['radial2'])
    def test_no_core_only_operator_replacement(self):self.assertFalse(audit.FLAGS['core_only_C16_operator_bound_valid_on_actual_activated_collar'])
    def test_full_pressure_memory_changes_third_and_fifth_target(self):self.assertTrue(self.r['full_pressure_memory_control']['different_third_and_fifth_rows'])
    def test_global_moment_half_factors(self):
        a=w.moment_bounds(2,3,5,7,(11,13));b=w.moment_bounds(2,3,5,7,(0,0))
        self.assertEqual(a[2]-b[2],F(11,2));self.assertEqual(a[4]-b[4],F(13,2))
    def test_first_moment_and_radial_weight(self):self.assertEqual(w.moment_bounds(2,3,5,7,(11,13))[:2],(6,12))
    def test_fourth_current_moment_has_both_linear_cross_terms(self):self.assertEqual(w.moment_bounds(2,3,5,7,(0,0))[3],144)
    def test_full_reciprocal_amplitude_cost_four(self):self.assertEqual(w.repair_bound((0,0,3,0,0),2,4,F(1,100))['transformed'][-1],3)
    def test_repair_uses_800_linear_inverse(self):
        r=w.repair_bound((1,2,3,4,5),2,3,F(1,1000));self.assertEqual(r['coefficient'],800*max(r['transformed']))
    def test_small_patch_slope_cost_not_erased(self):
        a=w.repair_bound((1,2,3,4,5),2,3,F(1,100));b=w.repair_bound((1,2,3,4,5),2,3,F(1,1000))
        self.assertEqual(b['transformed'][1],10*a['transformed'][1])
    def test_amplitude_eta_cost_not_frozen(self):self.assertEqual(w.repair_bound((0,0,1,0,0),1,2,F(1,100))['transformed'][-1],2)
    def test_exact_E_to_F_first_derivative(self):self.assertTrue(self.r['exact_radial_and_E_to_F_identities']['E_to_F_first'])
    def test_exact_E_to_F_second_derivative(self):self.assertTrue(self.r['exact_radial_and_E_to_F_identities']['E_to_F_second'])
    def test_logarithmic_second_derivative_keeps_minus_first(self):self.assertTrue(self.r['exact_radial_and_E_to_F_identities']['logarithmic_bump_second'])
    def test_cutoff_first_product_cost(self):
        c=w.collar_jets(1,1,F(1,2),1);r=w.extended_jets(c,0,1,3,5,2,1,(1,2,3));self.assertEqual(r['radial'],c['radial']+3*c['value'])
    def test_cutoff_second_product_cost(self):
        c=w.collar_jets(1,1,F(1,2),1);r=w.extended_jets(c,0,1,3,5,2,1,(1,2,3));self.assertEqual(r['radial2'],c['radial2']+6*c['radial']+5*c['value'])
    def test_patch_F_second_derivative_three_quarters_term(self):
        c=w.collar_jets(1,1,F(1,2),1);r=w.extended_jets(c,1,1,0,0,1,1,(1,0,0));self.assertEqual(r['radial2'],c['radial2']+F(3,4))
    def test_patch_lower_radius_enclosure_required(self):
        with self.assertRaises(ValueError):w.extended_jets(w.collar_jets(1,1,1,1),1,1,1,1,1,2,(1,1,1))
    def test_exact_theta_integral_measure(self):self.assertTrue(self.r['exact_stress_Jacobians']['theta_measure'])
    def test_exact_axial_integral_measure(self):self.assertTrue(self.r['exact_stress_Jacobians']['axial_measure'])
    def test_missing_theta_X_weight_detected(self):self.assertTrue(self.r['exact_stress_Jacobians']['missing_theta_X_weight_detected'])
    def test_missing_axial_sqrt_detected(self):self.assertTrue(self.r['exact_stress_Jacobians']['missing_axial_sqrt_detected'])
    def test_all_ten_exact_source_triangle_controls(self):self.assertEqual(len(self.r['source_triangle_controls']),10);self.assertTrue(all(all(x['checks'].values()) for x in self.r['source_triangle_controls']))
    def test_endpoint_angular_controls_present(self):self.assertEqual({r['eta'] for r in self.r['source_triangle_controls']},{'-1','-1/3','0','2/5','1'})
    def test_pressure_axial_derivative_cannot_be_omitted(self):self.assertTrue(all(r['checks']['pressure_enters_axial'] for r in self.r['source_triangle_controls']))
    def test_leading_axial_viscosity_charged_to_both_components(self):
        args=(F(1,200),1,1,2,dict(value=3,radial=4,radial2=5),7,11,13,17,19,23)
        a=w.source_bounds(*args,29,31);b=w.source_bounds(*args,0,0)
        self.assertEqual(a['theta']-b['theta'],29);self.assertEqual(a['axial']-b['axial'],31)
    def test_pressure_value_from_zero_axis_integral_not_just_local_derivative(self):
        r=w.source_bounds(F(1,200),1,1,2,dict(value=3,radial=4,radial2=5),7,11,13,17,19,23,29,31)
        self.assertEqual(r['pressure'],4*r['pressure_X']);self.assertEqual(r['pressure_axial'],5*r['pressure']/F(99,100))
    def test_independent_complete_first_series_comparison(self):
        r=self.r['independent_series_diagnostic'];self.assertEqual(r['case_count'],54);self.assertTrue(all(r['checks'].values()))
    def test_compact_stress_floor_never_zero(self):
        with self.assertRaises(ValueError):w.compact_gate(1,2,dict(value=1,radial=1),dict(theta=1,axial=1),0)
    def test_smaller_actual_stress_floor_increases_W_bound(self):
        args=(1,2,dict(value=1,radial=1),dict(theta=1,axial=1))
        self.assertEqual(w.compact_gate(*args,F(1,10))['A_W'],10*w.compact_gate(*args,1)['A_W'])
    def test_composed_Cstar_exactly_matches_previous_gate_formula(self):
        p=Path(__file__).resolve().parent.parent/'regional_weighted_gate/weighted_gate.py'
        s=importlib.util.spec_from_file_location('first_parent_gate',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
        j=dict(value=3,radial=4);r=w.compact_gate(1,2,j,dict(theta=5,axial=7),F(1,10))
        self.assertEqual(r['Cval'],m.value_constant(r['A_E'],r['A_EX'],r['A_UX'],r['A_W'],1,4)['Cval'])
    def test_invalid_collar_inputs_rejected(self):
        for args in ((0,1,1,1),(1,0,1,1),(1,1,0,1),(1,1,2,1),(1,1,1,0)):
            with self.assertRaises(ValueError):w.collar_jets(*args)
    def test_invalid_linear_normalizers_rejected(self):
        for r,e,lam in ((0,1,F(1,100)),(1,0,F(1,100)),(1,1,0),(1,1,F(1,2))):
            with self.assertRaises(ValueError):w.repair_bound((1,)*5,r,e,lam)
    def test_float_and_bool_are_not_exact_primitive_bounds(self):
        for x in (1.,True):
            with self.assertRaises(TypeError):w.collar_jets(1,x,1,1)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'b29ca142b37c4200d5a90b824cb2824f70cc0688');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_original_sources_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_unsupported_claims_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_primitive_quantitative_obligation_not_removable(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_boolean_integer_tamper_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_constant_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['rational_propagation_example']['gate']['Cval']='1'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
