"""Exact amplitude-join controls, frozen provenance and unsupported-scope guards."""
from copy import deepcopy
from fractions import Fraction as F
import importlib.util,json,tempfile,unittest
from pathlib import Path
import audit
import floor_bounds as w

class RegionalStressFloorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_actual_inner_offset_not_moved(self):self.assertEqual(self.r['inner_offset_and_hierarchy']['inner_start'],'Xa*exp(tc/16)')
    def test_actual_inner_flat_exponent(self):self.assertEqual(self.r['inner_offset_and_hierarchy']['flat_exponent'],25600)
    def test_polynomial_factor_two_not_dropped(self):self.assertIn('C^-25602',self.r['inner_offset_and_hierarchy']['inner_stress_lower'])
    def test_hierarchy_power_comparison(self):self.assertEqual(self.r['inner_offset_and_hierarchy']['exponent_difference'],15358)
    def test_inner_comparison_prefactor(self):self.assertTrue(self.r['inner_offset_and_hierarchy']['checks']['comparison_prefactor_correct'])
    def test_interval_begins_strictly_after_flat_zero(self):self.assertTrue(self.r['inner_offset_and_hierarchy']['checks']['fixed_interval_begins_strictly_after_flat_zero'])
    def test_log_hierarchy_uses_exact_integer_coefficient(self):self.assertTrue(self.r['inner_offset_and_hierarchy']['checks']['polynomial_log_comparison_coefficient'])
    def test_all_81_compact_controls_have_valid_inputs(self):
        rows=self.r['compact_coordinate_floor']['rows'];self.assertEqual(len(rows),81);self.assertTrue(all(r['checks']['input_bounds'] for r in rows))
    def test_all_rotated_norm_identities(self):self.assertTrue(all(r['checks']['exact_rotated_norm'] for r in self.r['compact_coordinate_floor']['rows']))
    def test_slope_ratio_budget_and_square_root(self):self.assertTrue(all(r['checks']['slope_ratio'] and r['checks']['square_root_denominator'] for r in self.r['compact_coordinate_floor']['rows']))
    def test_compact_D_margin_gives_amplitude_not_just_sign(self):self.assertTrue(all(r['checks']['compact_floor'] for r in self.r['compact_coordinate_floor']['rows']))
    def test_middle_floor_is_exactly_same_hierarchy_constant(self):
        for G in (F(1),F(2),F(19)):self.assertEqual(w.compact_lower(F(1,7),G**-4,4*G**3),w.common_floor(F(1,7),G))
    def test_stress_identity_matches_original_bridge(self):
        p=Path(__file__).resolve().parent.parent/'compact_pressure_bridge/bridge.py'
        s=importlib.util.spec_from_file_location('regional_floor_bridge',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
        for a,b,p1,p2 in ((3,0,4,0),(3,1,5,-1),(3,-1,5,1)):
            r=w.compact_coordinates(F(2,5),a,b,p1,p2);q=m.pulse(F(a),F(b),F(p1),F(p2),F(2,5))
            self.assertEqual(r['direct_stress2'],q['stress2']);self.assertEqual(r['rotated_stress2'],q['stress2_rotated'])
    def test_all_48_compensated_exterior_controls(self):
        rows=self.r['matched_compensated_exterior']['rows'];self.assertEqual(len(rows),48);self.assertTrue(all(all(r['checks'].values()) for r in rows))
    def test_exterior_shear_edit_is_charged(self):self.assertTrue(all(r['checks']['edited_shear_below_five'] for r in self.r['matched_compensated_exterior']['rows']))
    def test_exterior_finite_radius_not_sent_to_infinity(self):self.assertTrue(all(r['checks']['finite_radius_p_floor'] for r in self.r['matched_compensated_exterior']['rows']))
    def test_exterior_angular_stress_suffices(self):self.assertTrue(all(r['checks']['angular_stress_excess_exceeds_one'] for r in self.r['matched_compensated_exterior']['rows']))
    def test_same_pressure_and_functional_moment_prerequisite_explicit(self):self.assertIn('same forward stress by exact functional five-moment restoration',self.r['matched_compensated_exterior']['assumptions'])
    def test_all_64_first_bound_substitutions(self):
        rows=self.r['first_coefficient_substitution']['rows'];self.assertEqual(len(rows),64);self.assertTrue(all(all(r['checks'].values()) for r in rows))
    def test_reciprocal_stress_cost_not_lost(self):self.assertEqual(w.first_relative_bound(1,2,3,5,7),64*2**10*(F(9,2)*5+3*7))
    def test_previous_first_bound_formula_same_after_substitution(self):
        p=Path(__file__).resolve().parent.parent/'first_seminorm_bound/seminorms.py'
        s=importlib.util.spec_from_file_location('regional_floor_first_bound',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
        tau=w.common_floor(F(1,7),2);r=m.compact_gate(1,3,dict(value=5,radial=7),dict(theta=11,axial=13),tau)
        self.assertEqual(r['A_W'],w.first_relative_bound(F(1,7),2,3,11,13))
    def test_swirl_and_growing_shear_not_stress_lower_bound(self):self.assertTrue(self.r['scope_and_zero_edge_controls']['checks']['swirl_floor_alone_not_stress_floor'])
    def test_large_positive_stress_not_same_tilt_cone_certificate(self):self.assertGreater(F(self.r['scope_and_zero_edge_controls']['positive_large_stress_but_bad_cone_ratio2']),1)
    def test_flat_edge_not_in_floor_domain(self):self.assertTrue(self.r['scope_and_zero_edge_controls']['checks']['flat_endpoint_has_no_positive_unweighted_floor'])
    def test_logarithmic_underflow_control(self):self.assertTrue(self.r['scope_and_zero_edge_controls']['checks']['no_floating_underflow_interpreted_as_zero'])
    def test_log_floor_without_forming_G_to_tenth(self):self.assertEqual(w.log_floor(-100,10000)['logarithm_before_minus_log64'],-100100)
    def test_floor_scales_with_swirl_and_inverse_tenth_power(self):self.assertEqual(w.common_floor(2,2),F(1,512)*w.common_floor(1,1))
    def test_zero_source_has_zero_first_relative_bound(self):self.assertEqual(w.first_relative_bound(F(1,7),2,3,0,0),0)
    def test_invalid_or_zero_floor_inputs_rejected(self):
        for args in ((0,1),(1,0),(1,F(1,2)),(-1,2)):
            with self.assertRaises(ValueError):w.common_floor(*args)
    def test_invalid_compact_bounds_rejected(self):
        for args in ((0,1,1),(1,0,1),(1,1,F(1,2))):
            with self.assertRaises(ValueError):w.compact_lower(*args)
    def test_invalid_radius_and_source_bounds_rejected(self):
        for args in ((1,2,0,1,1),(1,2,1,-1,1)):
            with self.assertRaises(ValueError):w.first_relative_bound(*args)
    def test_float_and_bool_not_exact(self):
        for x in (1.,True):
            with self.assertRaises(TypeError):w.common_floor(x,1)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'5136b026a5626b21a283fba4f74a9ee424d66060');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record_reproduced(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_sources_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_no_scope_promotion(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_same_tilt_and_jet_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(1)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_boolean_integer_tamper_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_flat_exponent_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['inner_offset_and_hierarchy']['flat_exponent']=100
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_no_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
