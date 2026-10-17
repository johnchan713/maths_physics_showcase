"""Continuum inequalities, independent map controls and evidence integrity."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from matching_bounds import annular_review,target_review,repair_review,repair_budget_review,pressure_review,exp_bounds,vandermonde_rows
from matching_checks import polynomial_review,numerical_review,support_review
from audit import HERE,FLAGS,generate,validate,compare_records,protocol,parent_review,main
class MatchingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=generate()
    def test_01_parent(self):self.assertEqual(protocol()['parent_commit'],'655f9510308d39d001415d176bf259e085baf638')
    def test_02_protocol(self):self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
    def test_03_U_inverse(self):self.assertLess(max(map(F,annular_review()['U_inverse_rows'])),76)
    def test_04_E_inverse(self):self.assertLess(max(map(F,annular_review()['E_inverse_rows'])),610)
    def test_05_annular_quartic_jets(self):self.assertEqual(len(annular_review()['quadratic_normalized_eta_jet_upper']),5)
    def test_06_annular_q(self):self.assertLess(F(annular_review()['quadratic_algebra_upper']),50)
    def test_07_entry_logarithmic_absorption(self):self.assertTrue(target_review()['checks']['actual_current_entry_estimate_retained'])
    def test_08_target_all_rows(self):self.assertTrue(all(F(v)<F('.286') for v in target_review()['target_over_epsilon']))
    def test_09_target_keeps_pressure_row(self):self.assertEqual(F(target_review()['target_over_epsilon'][4]),F(1,4))
    def test_10_radius(self):self.assertEqual(F(target_review()['radius_over_epsilon']),F('348.92'))
    def test_11_contraction(self):self.assertLess(F(target_review()['contraction_parameter']),F(1,10**8))
    def test_12_same_C2_root(self):self.assertTrue(target_review()['checks']['same_old_C2_root'])
    def test_13_same_C1_root(self):self.assertTrue(target_review()['checks']['same_old_C1_root'])
    def test_14_repair_U_inverse(self):self.assertLess(F(repair_review()['U_inverse_upper']),2)
    def test_15_repair_E_inverse(self):self.assertLess(F(repair_review()['E_inverse_upper']),800)
    def test_16_repair_q(self):self.assertLess(F(repair_review()['quadratic_algebra_upper']),25000)
    def test_17_lambda_loss(self):self.assertIn('lambda*N',repair_budget_review()['target'])
    def test_18_same_frequency(self):self.assertIn('floor(H^32)',repair_budget_review()['frequency'])
    def test_19_full_sum_root(self):self.assertIn('1600*D/N',repair_budget_review()['root'])
    def test_20_logarithmic_contraction(self):self.assertLess(8*800**2*25000,2**37);self.assertLess(F(2)**-219,1)
    def test_21_original_e_normalization(self):self.assertTrue(repair_review()['checks']['constant_eta_B_and_Q_requires_original_e'])
    def test_22_Catalan_independent(self):self.assertTrue(polynomial_review()['checks']['independent_Catalan_and_implicit_series_agree'])
    def test_23_variable_Q_residual(self):self.assertTrue(polynomial_review()['checks']['exact_variable_Q_residual_vanishes'])
    def test_24_constant_Q_negative_control(self):self.assertTrue(polynomial_review()['checks']['constant_Q_jet_rule_fails'])
    def test_25_jet_majorant(self):self.assertTrue(polynomial_review()['checks']['triple_majorant_controls_all_eight_jets'])
    def test_26_factorial_jet_normalization(self):self.assertTrue(polynomial_review()['checks']['raw_second_derivative_needs_factor_two'])
    def test_27_energy_cancellation(self):self.assertTrue(polynomial_review()['checks']['restoration_energy_identity_exact'])
    def test_28_missing_pressure_negative_control(self):self.assertTrue(polynomial_review()['checks']['pressure_value_zero_does_not_remove_gradient'])
    def test_29_actual_annular_map(self):r=numerical_review();self.assertLess(r['annular_inverse'],610);self.assertLess(r['annular_residual'],1e-60)
    def test_30_actual_implicit_derivative(self):self.assertLess(numerical_review()['annular_eta_error'],1e-35)
    def test_31_actual_repair_maps(self):self.assertEqual(len(numerical_review()['repair_rows']),2);self.assertTrue(all(r['root_residual']<1e-60 and r['inverse']<800 for r in numerical_review()['repair_rows']))
    def test_32_ordinary_integrands(self):self.assertTrue(all(r['ordinary_control_error']<1e-60 for r in numerical_review()['repair_rows']))
    def test_33_lost_input_negative_control(self):self.assertTrue(numerical_review()['checks']['rounded_input_difference_can_be_lost'])
    def test_34_zero_lambda_rejected(self):self.assertTrue(numerical_review()['checks']['zero_lambda_rejected_despite_weight_extension'])
    def test_35_support_and_norm_controls(self):self.assertTrue(all(support_review()['checks'].values()))
    def test_36_all_pressure_rows_and_remaining_scope(self):self.assertEqual(pressure_review()['rows'],['M','J','I','S','Cp']);self.assertFalse(FLAGS['complete_core_global_backward_pressure_identification_discharged'])
    def test_37_bindings_and_frozen_record(self):
        self.assertEqual(len(self.r['inherited_provenance']),1069);self.assertEqual(len(self.r['source_hashes']),18)
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.r,f);validate(json.loads(json.dumps(self.r)))
    def test_38_scope_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.r);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_39_history_and_review_mutation_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        r=generate();r['reviews']['actual_annular_inverse_and_quartic_bound']['checks']['inverse_below_610']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1]);self.assertTrue(annular_review()['checks']['inverse_below_610'])
    def test_40_protocol_gate_source_scope_and_overwrite_guards(self):
        for k in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.r)
            if isinstance(r[k],dict):r[k].pop(next(iter(r[k])))
            else:r[k].pop()
            with self.assertRaises(ValueError):validate(r)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
if __name__=='__main__':unittest.main()
