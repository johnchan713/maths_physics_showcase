"""Actual native crossing, inverse pullback and necessary-only scope guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import pullback as p

class PullbackTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_actual_a_two(self):self.assertTrue(self.r['symbolic']['checks']['actual_a_is_two'])
    def test_exact_actual_axial_shear(self):self.assertTrue(self.r['symbolic']['checks']['actual_axial_shear_identity'])
    def test_taper_center_derivative(self):self.assertTrue(self.r['symbolic']['checks']['actual_taper_center_derivative'])
    def test_no_zero_pressure_replacement(self):self.assertIs(self.r['symbolic']['pressure_coordinate_set_to_zero'],False)
    def test_exact_four_inverse_derivatives(self):self.assertTrue(self.r['symbolic']['checks']['inverse_four_derivatives_exact'])
    def test_independent_implicit_inverse_four_orders(self):self.assertTrue(self.r['symbolic']['checks']['implicit_normalized_inverse_all_four_orders'])
    def test_normalized_inverse_constants(self):self.assertEqual(self.r['symbolic']['normalized_inverse_bounds'],[1,3,16,105])
    def test_exact_crossing_derivative(self):self.assertTrue(self.r['symbolic']['checks']['exact_angular_crossing_derivative'])
    def test_exact_entire_window_pullback(self):self.assertTrue(self.r['symbolic']['checks']['exact_window_pullback'])
    def test_exact_target_squared_identity(self):self.assertTrue(self.r['symbolic']['checks']['actual_target_squared_pullback'])
    def test_q_norm_not_raw_derivative_max(self):self.assertEqual(self.r['rational']['q_C4_over_r_upper'],'221/256');self.assertTrue(self.r['rational']['checks']['q_C4_below_r'])
    def test_full_inverse_norm_below_one(self):self.assertTrue(self.r['rational']['checks']['inverse_norm_below_one_at_r_bound'])
    def test_inverse_product_whole_norm126(self):self.assertTrue(self.r['rational']['checks']['inverse_whole_C4_cost_fits126'])
    def test_reciprocal_sqrt_constant(self):self.assertTrue(self.r['rational']['checks']['reciprocal_sqrt_C4_below_two'])
    def test_raw_shape_normalization_factor24(self):self.assertTrue(self.r['rational']['checks']['raw_fourth_factorial_is24'])
    def test_final_shape_extraction_factor192(self):self.assertTrue(self.r['rational']['checks']['shape_extraction_constant_192'])
    def test_outward_parameter_two_precisions(self):self.assertEqual([r['interval_digits'] for r in self.r['parameters']],[80,110]);self.assertTrue(all(all(r['checks'].values()) for r in self.r['parameters']))
    def test_slice_before_final_flat_units(self):self.assertTrue(all(r['checks']['slice_before_final_eleven_flat_units'] for r in self.r['parameters']))
    def test_actual_P1_below_A(self):self.assertTrue(all(r['checks']['P1_upper_below_A'] for r in self.r['parameters']))
    def test_actual_beta_below_one(self):self.assertTrue(all(r['checks']['beta_upper_below_one'] for r in self.r['parameters']))
    def test_no_enormous_actual_parameter_materialization(self):self.assertTrue(all(r['actual_A_beta_delta_G_H_N_or_radius_materialized'] is False for r in self.r['parameters']))
    def test_six_independent_inverse_and_fourth_diagnostics(self):self.assertEqual(sum(len(r['rows']) for r in self.r['numerical']),6);self.assertTrue(all(all(all(v['checks'].values()) for v in r['rows']) for r in self.r['numerical']))
    def test_both_diagnostic_scales(self):self.assertEqual([r['delta'] for r in self.r['numerical']],['1e-18','1e-24'])
    def test_manufactured_not_actual_or_interval(self):self.assertTrue(all(r['manufactured_floating_not_interval_or_actual_parameter_evaluation'] is True for r in self.r['numerical']))
    def test_frozen_tolerance_not_roundoff_string(self):self.assertTrue(all(r['relative_error_limit']=='1e-55' and 'relative_error' not in r for r in self.r['numerical']))
    def test_necessary_gate_can_reject_large_auxiliary_shape(self):self.assertFalse(p.necessary_gate(2**16,2**513)['necessary_condition_not_violated'])
    def test_passing_necessary_gate_does_not_accept_S(self):
        r=p.necessary_gate(2**16,0);self.assertTrue(r['necessary_condition_not_violated']);self.assertIs(r['actual_S_accepted'],False);self.assertIs(r['actual_shape_enclosed'],False)
    def test_gate_never_disproves_actual_H(self):self.assertIs(p.necessary_gate(2**16,2**513)['actual_H_disproved'],False)
    def test_invalid_gate_parameters(self):
        for A,Z in ((1,0),(True,0),(2**16,-1),(2**16,True),(2**16,0.5)):
            with self.assertRaises(ValueError):p.necessary_gate(A,Z)
    def test_diagnostic_domain_guard(self):
        for d,b in (('0.01','0.1'),('1','0.1'),('-1','0.1')):
            with self.assertRaises(ValueError):p.manufactured_review(d,b)
    def test_interval_precision_guard(self):
        with self.assertRaises(ValueError):p.parameter_review(40)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'8c34b682cd1d1ee593f190ae46efabe65a647f1f');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_fixed_historical_manifest_not_future_glob(self):self.assertEqual(len(audit.HISTORICAL_RECORDS),129);self.assertEqual(len(set(audit.HISTORICAL_RECORDS)),129);self.assertNotIn(audit.HERE.name,audit.HISTORICAL_RECORDS)
    def test_all_historical_bindings(self):self.assertIs(self.record['inherited_provenance']['all_bindings_verified'],True);self.assertGreater(self.record['inherited_provenance']['binding_count'],2228)
    def test_frozen_exact_reproduction(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_cannot_replace_false(self):
        r=deepcopy(self.record);r['flags']['actual_selected_Z_enclosed_or_S_condition_verified']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_necessary_sufficient_and_memory_protocol_cannot_be_erased(self):
        for k in ('native','geometry','pullback','bound','necessary','sufficient_retained','memory','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_actual_Z_obligation_not_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)]);audit.validate(json.loads(out.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(out)])

if __name__=='__main__':unittest.main()
