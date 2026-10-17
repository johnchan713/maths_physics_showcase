"""Exact frame/floor controls and retained scientific limits."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,main
from frame_floors import squared_floors,leading_identities,identity_review,global_review,logarithmic_review
from frame_checks import exact_frame,frame_review,nearest_nonzero,carrier_review,pressure_and_domain_review,matmul


class LeadingFrameTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()
    def test_01_protocol(self):self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
    def test_02_current_parent(self):self.assertEqual(protocol()['parent_commit'],'889f2e5ba0837a3ea90abdacd7ae37d92f366b3d')
    def test_03_global_closed_regions(self):self.assertEqual(len(global_review()['regions']),5);self.assertTrue(all(global_review()['checks'].values()))
    def test_04_unambiguous_half_gaps(self):self.assertIn('1/(2G^3)',protocol()['regional']);self.assertIn('1/(2G^4)',protocol()['regional'])
    def test_05_exact_shear_identities(self):self.assertTrue(all(identity_review()['checks'].values()));self.assertEqual(len(identity_review()['rows']),24)
    def test_06_cancellation_in_c0(self):self.assertEqual(leading_identities(2,3,4)['c2'],F(19,6))
    def test_07_growth_zero_at_v_two(self):self.assertEqual(leading_identities(1,2,0)['lambda2'],0)
    def test_08_a_positive_not_growth(self):self.assertLess(leading_identities(1,1,0)['lambda2'],0)
    def test_09_zero_swirl_rejected(self):
        with self.assertRaises(ValueError):leading_identities(0,3,0)
    def test_10_invalid_squared_floors(self):
        for p in ((0,2,1,4,1,F(1,2)),(1,2,0,4,1,F(1,2)),(1,2,1,4,1,0),(2,1,1,4,1,F(1,2)),(1,2,1,1,0,F(1,2))):
            with self.assertRaises(ValueError):squared_floors(*p)
    def test_11_frame_inverse_is_inverse_gap(self):self.assertEqual(squared_floors(1,2,1,4,2,F(1,8))['left_inverse_max2'],8)
    def test_12_strong_gap_inverse_saturates(self):self.assertEqual(squared_floors(1,2,1,8,2,4)['left_inverse_max2'],F(1,2))
    def test_13_log_floor_survives_underflow(self):self.assertTrue(logarithmic_review()['checks']['floating_underflow_not_a_zero_analytical_floor'])
    def test_14_nonconstant_logarithmic_profile(self):r=logarithmic_review();self.assertLess(r['log_floor'],r['log_actual'])
    def test_15_open_positivity_not_uniform_floor(self):self.assertTrue(logarithmic_review()['checks']['pointwise_open_positivity_does_not_give_closed_uniform_floor'])
    def test_16_actual_ambient_left_inverse(self):self.assertTrue(frame_review()['checks']['specified_ambient_left_inverse_not_silently_replaced'])
    def test_17_frame_transversality_and_products(self):self.assertEqual(frame_review()['exact_frame_cases'],12);self.assertTrue(all(frame_review()['checks'].values()))
    def test_18_degenerate_frame_rejected(self):
        for t,s,c,w in (((1,1),0,1,1),((1,0),0,0,1),((1,0),1,1,1)):
            with self.assertRaises(ValueError):exact_frame(t,s,c,w)
    def test_19_original_frame_left_product(self):
        r=exact_frame((F(3,5),F(4,5)),F(3,4),F(-1,2),F(5,4));self.assertEqual(matmul(r['left'],r['B']),[[1,0],[0,1]])
    def test_20_round_zero_nonzero_and_full_error(self):self.assertNotEqual(nearest_nonzero(0,10),0);self.assertEqual(abs(nearest_nonzero(0,10)),F(1,10))
    def test_21_rounding_half_bound_fails(self):self.assertTrue(carrier_review()['checks']['half_over_k_bound_fails_at_zero'])
    def test_22_carrier_both_signs(self):self.assertEqual(len(carrier_review()['rows']),8);self.assertTrue(all(carrier_review()['checks'].values()))
    def test_23_disappearing_angular_component(self):self.assertTrue(carrier_review()['checks']['angular_component_can_tend_to_zero_with_nonzero_integer_frequency'])
    def test_24_axial_floor_not_moving_normal(self):self.assertEqual(carrier_review()['axial_normalized_upper'],'-1/8');self.assertTrue(carrier_review()['checks']['moving_axial_normal_not_identified_with_constant_pz'])
    def test_25_invalid_rounding_integer(self):
        for k in (0,-1,1.5):
            with self.assertRaises(ValueError):nearest_nonzero(1,k)
    def test_26_pressure_transversality_is_not_zero_pressure(self):self.assertEqual(pressure_and_domain_review()['normal_pressure_numerator_without_K'],'3');self.assertTrue(all(pressure_and_domain_review()['checks'].values()))
    def test_27_retained_tilt_and_supports(self):self.assertIn('original tilt',protocol()['fixed_choices']);self.assertIn('no new support',protocol()['support'])
    def test_28_no_unconditional_promotion(self):self.assertFalse(FLAGS['unconditional_pulse_family_or_blowup_proved']);self.assertFalse(FLAGS['every_generated_complete_covariance_mean_and_residual_input_closed'])
    def test_29_source_record_and_parent_history(self):
        old=json.loads((HERE/'evidence.json').read_text());validate(old);compare_records(self.record,old)
        self.assertEqual(len(self.record['source_hashes']),19);self.assertGreater(len(self.record['inherited_provenance']),3100)
    def test_30_reject_false_scope_flags(self):
        for key,value in FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):validate(r)
    def test_31_reject_missing_history_and_limits(self):
        for key in ('source_hashes','inherited_provenance','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[key],dict):r[key].pop(next(iter(r[key])))
            else:r[key].pop()
            with self.assertRaises(ValueError):validate(r)
    def test_32_reject_changed_floor_formula(self):
        for key in ('spectral','rounding','frame'):
            r=deepcopy(self.record);r['protocol'][key]='zero'
            with self.assertRaises(ValueError):validate(r)
    def test_33_reject_failed_control(self):
        r=deepcopy(self.record);r['reviews']['original_moving_frame_and_ambient_left_inverse']['checks']['original_left_inverse_identity']=False
        with self.assertRaises(ValueError):validate(r)
    def test_34_fresh_import_preserves_paths(self):
        code='import sys;sys.path.insert(0,'+repr(str(HERE))+');before=sys.path[:];import audit;assert audit.HERE==__import__("pathlib").Path('+repr(str(HERE))+');assert before==sys.path;assert "frame_floors" in sys.modules'
        subprocess.run([sys.executable,'-c',code],check=True,capture_output=True,text=True)
    def test_35_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_36_roundtrip_and_all_gates(self):self.assertTrue(all(self.record['gates'].values()));validate(json.loads(json.dumps(self.record)))
if __name__=='__main__':unittest.main()
