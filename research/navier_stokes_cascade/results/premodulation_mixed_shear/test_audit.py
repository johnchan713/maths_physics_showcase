"""Native mixed-input estimates, variance-shape obstruction and scope guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import mixed_bounds as m

class MixedShearTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_annular_relative_factor_positive(self):self.assertLess(m.annular_bounds()['relative_C4'],Q(1,2))
    def test_annular_actual_Ey_C4_below_P(self):self.assertLess(m.annular_bounds()['Ey_C4_over_P'],1)
    def test_annular_actual_Uy_C4_below_one(self):self.assertLess(m.annular_bounds()['Uy_C4'],1)
    def test_annular_physical_a_C4_below_one(self):self.assertLess(m.annular_bounds()['physical_a_C4'],1)
    def test_earlier_width_point08_not_leading_point1(self):self.assertEqual(m.annular_bounds()['leading_bump_width'],Q(2,25));self.assertNotEqual(m.annular_bounds()['leading_bump_width'],Q(1,10))
    def test_all_three_E_and_two_U_coefficients_charged(self):
        q=m.annular_bounds();d=Q(572,10**16);self.assertEqual(q['relative_y_C4'],3*d*250000);self.assertEqual(q['Uy_C4'],2*d*250000)
    def test_relative_inverse_rejects_large_unvalidated_root(self):
        with self.assertRaises(ValueError):m.annular_bounds(d=1)
    def test_negative_annular_bounds_rejected(self):
        with self.assertRaises(ValueError):m.annular_bounds(b1=-1)
    def test_exact_annular_inputs_required(self):
        with self.assertRaises(TypeError):m.annular_bounds(d=5.72e-14)
    def test_shape_slope_below_one(self):self.assertLess(m.regional_bounds(682)['shape_log_slope_C4'],1)
    def test_outer_interpolation_derivative_is_charged(self):self.assertEqual(m.regional_bounds(682)['outer_log_slope_C4'],1+Q(64,200)*683)
    def test_cover_sums_regions_not_maximum_norm(self):
        q=m.regional_bounds(682);self.assertGreater(q['log_slope_C4_cover_sum'],max(682,1000));self.assertLess(q['log_slope_C4_cover_sum'],q['log_slope_C4_majorant'])
    def test_all_axial_regions_charged(self):self.assertEqual(m.regional_bounds(1000)['Uy_C4_cover_sum'],1010)
    def test_regional_minimum_input_checked(self):
        with self.assertRaises(ValueError):m.regional_bounds(681)
    def test_two_interval_precisions(self):self.assertEqual([v['interval_digits'] for v in self.r['parameters']],[80,110])
    def test_actual_reference_p1_C4_comparison(self):self.assertTrue(all(v['checks']['reference_p1_C4_below_Bold_over100'] for v in self.r['parameters']))
    def test_actual_prescribed_axial_bound_comparison(self):self.assertTrue(all(v['checks']['actual_prescribed_Uy_C4_below_Bold'] for v in self.r['parameters']))
    def test_global_native_a_and_Uy_below_C(self):self.assertTrue(all(v['checks']['global_a_C4_below_C'] and v['checks']['global_Uy_C4_below_C'] for v in self.r['parameters']))
    def test_shape_precedes_actual_restoration(self):self.assertTrue(all(v['checks']['shape_exit_precedes_restoration'] for v in self.r['parameters']))
    def test_large_frequency_and_profile_not_materialized(self):self.assertTrue(all(v['actual_enormous_C_G_H_N_or_modulated_profile_materialized'] is False for v in self.r['parameters']))
    def test_actual_prescribed_E_identity(self):self.assertTrue(self.r['symbolic']['checks']['actual_prescribed_E_radial_identity'])
    def test_actual_prescribed_axial_identity(self):self.assertTrue(self.r['symbolic']['checks']['actual_prescribed_U_radial_identity'])
    def test_mixed_eta4_E_needs_no_pure_fifth(self):self.assertTrue(self.r['symbolic']['checks']['actual_mixed_E_eta4_uses_only_slope_eta4'])
    def test_mixed_eta4_U_needs_no_pure_fifth(self):self.assertTrue(self.r['symbolic']['checks']['actual_mixed_U_eta4_uses_only_slope_eta4'])
    def test_reference_shear_not_copied_after_flattening(self):self.assertTrue(self.r['symbolic']['checks']['flattened_reference_not_substituted_for_actual_prescribed_shear'])
    def test_relative_bump_a_identity(self):self.assertTrue(self.r['symbolic']['checks']['relative_annular_physical_a_identity'])
    def test_relative_bump_mixed_E_identity(self):self.assertTrue(self.r['symbolic']['checks']['relative_annular_E_mixed_identity'])
    def test_additive_swirl_mistake_detected(self):self.assertTrue(self.r['symbolic']['checks']['treating_relative_E_bump_as_additive_detected'])
    def test_native_loop_algebra_exponents(self):
        q=m.input_exponents();self.assertEqual((q['a_C4'],q['inverse_a_C4'],q['bs_C4'],q['ts_C4'],q['vs_C4']),(1,10,1286,1296,2583))
    def test_variance_target_full_axial_square_power(self):self.assertEqual(m.input_exponents()['sqrt_target_C4_A_power'],8*2583+9)
    def test_target_G_absorption_is_conditional_on_shape(self):self.assertTrue(all(v['checks']['same_G_absorbs_fourth_cutoff_template_if_shape_Z_le_1e100'] for v in self.r['parameters']))
    def test_sqrt_fourth_coefficients_retained(self):self.assertTrue(self.r['symbolic']['checks']['sqrt_fourth_raw_coefficient_is_minus15_over16']);self.assertTrue(self.r['symbolic']['checks']['inverse_sqrt_fourth_raw_coefficient_is105_over16'])
    def test_compressed_cutoff_third_center_for_fourth_MVT(self):self.assertEqual(self.r['symbolic']['canonical_step_third_at_center'],'-3712')
    def test_compressed_fourth_lower_bound_scales_n4(self):
        for n in (1,2,8,64):self.assertEqual(Q(3712*n**3)/Q(1,2*n),7424*n**4)
    def test_auxiliary_small_shape_gate_accepts_not_actual(self):
        r=m.cutoff_shape_gate(2**16,10**100);self.assertTrue(r['accepts']);self.assertFalse(r['actual_cutoff_budget_supplied']);self.assertTrue(r['auxiliary_only'])
    def test_shape_gate_can_reject(self):self.assertFalse(m.cutoff_shape_gate(2**16,1,logZ=str(2**513))['accepts'])
    def test_gate_domain_checked(self):
        for A,Z in ((2,1),(2**16,0),(True,1),(2**16,1.)):
            with self.assertRaises(ValueError):m.cutoff_shape_gate(A,Z)
    def test_manuscript_SHA_retained(self):self.assertEqual(audit.protocol()['manuscript_sha256'],'0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f')
    def test_parent_protocol_and_frozen_record(self):self.assertEqual(audit.protocol()['parent_commit'],'04a9058da056821a2f4833c10550ec958a9424c1');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_source_bindings(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_unproved_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_not_boolean_false(self):
        r=deepcopy(self.record);r['flags']['independent_foundation_review_or_formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_cutoff_correction_native_scope_and_target_not_erased(self):
        for k in ('source_correction','native','annular','cover','shape','target','gate','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_variance_shape_obligation_not_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_evidence_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
