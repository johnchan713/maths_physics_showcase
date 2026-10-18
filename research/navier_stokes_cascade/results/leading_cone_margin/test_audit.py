"""Exact direction, existing-tilt and old-gap-map controls, not actual profile samples."""
from copy import deepcopy
from fractions import Fraction as F
import importlib.util,json,tempfile,unittest
from pathlib import Path
import audit
import cone_bounds as w

class LeadingConeMarginTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_compact_orientation_constant(self):self.assertEqual(w.margins(2)['omega'],F(1,640*2**13))
    def test_compact_raw_deficit_constant(self):self.assertEqual(w.margins(2)['epsilon'],F(1,2**40))
    def test_all_81_compact_controls(self):
        rows=self.r['compact_direction_margins']['rows'];self.assertEqual(len(rows),81);self.assertTrue(all(all(r['checks'].values()) for r in rows))
    def test_full_fourth_gap_normalization(self):self.assertTrue(all(r['checks']['exact_raw_deficit_identity'] for r in self.r['compact_direction_margins']['rows']))
    def test_orientation_is_unit_direction_not_amplitude(self):
        r=w.coordinates(4,0,F(24,5),F(3,5));self.assertEqual(r['orientation2'],F(16,25));self.assertEqual(r['raw_ratio2'],F(9,16))
    def test_original_bridge_direction_identity(self):
        p=Path(__file__).resolve().parent.parent/'compact_pressure_bridge/bridge.py'
        s=importlib.util.spec_from_file_location('leading_cone_old_bridge',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
        r=w.coordinates(4,0,F(24,5),F(3,5));q=m.pulse(F(4),F(0),F(24,5),F(3,5),F(1))
        self.assertEqual(r['orientation2'],q['orientation2']);self.assertEqual(r['raw_ratio2'],q['ratio2'])
    def test_all_72_finite_p_exterior_controls(self):
        rows=self.r['compensated_exterior_direction_margins']['rows'];self.assertEqual(len(rows),72);self.assertTrue(all(all(r['checks'].values()) for r in rows))
    def test_exterior_finite_p_term_not_dropped(self):self.assertTrue(self.r['compensated_exterior_direction_margins']['checks']['exact_finite_p_lower'])
    def test_exterior_w_and_b_over_a_costs(self):self.assertTrue(self.r['compensated_exterior_direction_margins']['checks']['actual_w_and_shear_ratio_budget'])
    def test_exterior_raw_deficit_denominator(self):self.assertEqual(F(4,5)/(2*14000**2),F(1,490000000))
    def test_inner_orientation_hierarchy_exponent(self):self.assertTrue(self.r['closed_annulus_direction_join']['checks']['inner_orientation_power_dominance'])
    def test_global_constants_below_exterior_bounds(self):
        r=w.margins(2);self.assertLess(r['omega'],F(1,64000));self.assertLess(r['epsilon'],F(1,490000000))
    def test_terminal_bounds_extend_to_limiting_direction(self):self.assertTrue(self.r['closed_annulus_direction_join']['checks']['all_terminal_direction_bounds'])
    def test_fixed_tilt_formula_has_u_inverse_squared(self):
        r=w.tilt_gate(F(1,4),4);self.assertEqual(r['squared_ratio_upper'],F(15,16));self.assertEqual(r['delta'],F(1,16))
    def test_exact_sufficient_threshold(self):self.assertEqual(w.tilt_gate(F(1,4),4)['sharp_sufficient_u2_threshold'],3)
    def test_equality_rejected_for_strict_margin(self):self.assertFalse(w.tilt_gate(F(1,4),3)['accepted'])
    def test_below_threshold_returns_uncertified_not_negative_gamma(self):
        r=w.tilt_gate(F(1,4),1);self.assertFalse(r['accepted']);self.assertIsNone(r['gamma'])
    def test_safe_sqrt_margin(self):
        r=w.tilt_gate(F(1,4),4);self.assertEqual(r['gamma'],F(1,32));self.assertGreaterEqual((1-r['gamma'])**2,r['squared_ratio_upper'])
    def test_padded_tilt_gate(self):
        r=w.tilt_gate(F(1,4),8);self.assertTrue(r['padded_accepted']);self.assertEqual(r['padded_gamma'],F(1,16))
    def test_all_15_exact_tilt_controls(self):
        rows=self.r['existing_tilt_gate']['rows'];self.assertEqual(len(rows),15);self.assertTrue(all(all(r['checks'].values()) for r in rows))
    def test_failure_of_coarse_gate_not_failure_of_actual_cone(self):self.assertTrue(self.r['existing_tilt_gate']['checks']['coarse_gate_failure_is_not_actual_cone_failure'])
    def test_sharp_boundary_control(self):self.assertTrue(self.r['existing_tilt_gate']['checks']['equality_cannot_supply_strict_margin'])
    def test_boundaries_on_both_sides(self):self.assertTrue(self.r['existing_tilt_gate']['checks']['both_sides_of_sharp_boundary'])
    def test_chart_positive_scaling_keeps_direction(self):
        H,J,c2=F(4,5),F(3,5),F(1)
        for scale in (F(1,5),F(1),F(5)):self.assertEqual(c2*(scale*J)**2/(scale*H)**2,F(9,16))
    def test_chart_floor_and_upper_factors(self):
        r=w.chart_inputs(F(1,10),2,F(1,4),4,1,F(1,10));self.assertEqual(r['f'],F(1,50));self.assertEqual(r['Fbar'],10)
    def test_chart_squared_discriminant_factor(self):self.assertEqual(w.chart_inputs(F(1,10),2,F(1,4),4,1,F(1,10))['d'],F(1,100000))
    def test_parent_gap_radius_positive(self):self.assertTrue(self.r['scaled_original_gap_interface']['checks']['unchanged_parent_radius_is_positive'])
    def test_original_same_u_is_not_rechosen(self):self.assertTrue(audit.FLAGS['existing_tilt_frame_carriers_roots_frequency_width_and_supports_unchanged']);self.assertFalse(audit.FLAGS['actual_existing_tilt_satisfies_new_sufficient_inequality_verified'])
    def test_invalid_hierarchy_range(self):
        with self.assertRaises(ValueError):w.margins(1)
    def test_invalid_tilt_inputs(self):
        for e,U in ((0,1),(1,1),(F(1,4),0),(F(1,4),-1)):
            with self.assertRaises(ValueError):w.tilt_gate(e,U)
    def test_invalid_chart_envelopes(self):
        with self.assertRaises(ValueError):w.chart_inputs(0,2,1,4,1,1)
    def test_invalid_unoriented_coordinates(self):
        with self.assertRaises(ValueError):w.coordinates(3,0,3,0)
    def test_exact_inputs_reject_floats_and_bools(self):
        for x in (2.,True):
            with self.assertRaises(TypeError):w.margins(x)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'1f031a3c84cfafc8461bdb531d05adde904581b9');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_sources_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_no_scope_promotion(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_same_tilt_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(1)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_integer_boolean_scope_mutation_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_margin_constant_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['closed_annulus_direction_join']['epsilon']='1/2'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
