"""Exact cone/root controls, Haar diagnostics and scientific scope guards."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,parent_review,main
from covariance_bounds import cone_bounds,matrix,determinant,solve,freezing_margin,length_requirement,root_interval,corner_review,scale_review
from covariance_checks import reciprocal,square_root,polynomial,mul,jet_review,gaussian_review,edge_review


class PrimaryCovarianceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()
    def test_01_protocol(self):self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
    def test_02_current_parent(self):self.assertEqual(protocol()['parent_commit'],'b17351e44a98b0ca8ef02a936ea7a0d0e94a4fc4')
    def test_03_all_corner_and_target_cases(self):self.assertEqual(corner_review()['corner_target_cases'],384);self.assertTrue(corner_review()['checks']['all_rational_corner_targets_satisfy_derived_bounds'])
    def test_04_safe_normalized_inverse(self):self.assertEqual(cone_bounds(1,F(1,4))['inverse_upper'],F(20,9));self.assertLess(F(20,9),3)
    def test_05_safe_determinant(self):self.assertEqual(cone_bounds(1,F(1,4))['det_lower'],F(9,8))
    def test_06_safe_positive_roots(self):
        for g in (F(1,100),F(1,2),F(1)):
            b=cone_bounds(g,g/4);self.assertGreater(b['root_lower'],g/8);self.assertLess(b['root_upper'],3)
    def test_07_sharp_failure_despite_invertibility(self):
        g=F(1,2);critical=g/(2-g)
        for eta in (critical,critical+F(1,100)):
            a=matrix((0,eta,0,-eta));self.assertEqual(determinant(a),2);self.assertLessEqual(solve(a,(1,1-g))[0],0)
    def test_08_invalid_cone_inputs(self):
        for g,e in ((0,0),(2,0),(1,-1),(1,1)):
            with self.assertRaises(ValueError):cone_bounds(g,e)
    def test_09_exact_Cramer_solution(self):self.assertEqual(solve(matrix((0,0,0,0)),(F(3),F(1))),(1,2))
    def test_10_singular_solve_rejected(self):
        with self.assertRaises(ValueError):solve(((1,1),(1,1)),(1,1))
    def test_11_freezing_safe_budget(self):
        for g in (F(1,100),F(1,2),F(1)):self.assertGreaterEqual(freezing_margin(g,g/8),g/2)
    def test_12_exact_half_margin_threshold(self):
        for g in (F(1,10),F(1,2),F(1)):self.assertEqual(freezing_margin(g,g/(4-g)),g/2)
    def test_13_retained_length_budget_dominates(self):self.assertTrue(scale_review()['checks']['retained_pulse_length_budget_dominates_new_explicit_condition'])
    def test_14_invalid_length_or_freezing_inputs(self):
        with self.assertRaises(ValueError):length_requirement(1,0,1,F(1,2))
        with self.assertRaises(ValueError):freezing_margin(1,1)
    def test_15_root_bounds_keep_Haar_scale(self):
        lo,hi=root_interval(F(1,2),4,F(1,100),F(1,50));self.assertEqual(lo,F(25,2));self.assertEqual(hi,1200)
    def test_16_exact_Taylor_inverse(self):self.assertTrue(jet_review()['checks']['independent_Cramer_Taylor_inverse_exact'])
    def test_17_root_Taylor_jets(self):
        r=jet_review();self.assertEqual(r['normalized_Taylor_order'],6);self.assertTrue(r['checks']['differentiated_positive_square_roots_recovered_exactly'])
    def test_18_signed_map_uses_original_positive_root(self):
        r=jet_review()['checks'];self.assertTrue(r['original_root_signed_linear_covariance_identity']);self.assertTrue(r['updated_square_root_not_required_for_signed_target'])
    def test_19_signed_square_not_erased(self):self.assertTrue(jet_review()['checks']['signed_square_is_nonzero_and_not_erased'])
    def test_20_reciprocal_and_square_root_jet_controls(self):
        a=polynomial((2,1));self.assertEqual(mul(a,reciprocal(a)),polynomial((1,)));self.assertEqual(square_root(mul(a,a),2),a)
        with self.assertRaises(ValueError):square_root(polynomial((2,)),1)
    def test_21_all_Haar_mass_controls(self):self.assertEqual(len(gaussian_review()['rows']),4);self.assertTrue(gaussian_review()['checks']['actual_Haar_area_bounds_not_radius_free'])
    def test_22_gaussian_direction_budget(self):self.assertTrue(gaussian_review()['checks']['Gaussian_column_errors_fit_safe_cone_budget'])
    def test_23_homogeneous_Riccati_barriers(self):self.assertTrue(gaussian_review()['checks']['homogeneous_Riccati_radial_barriers'])
    def test_24_actual_outer_weighted_pole_not_suppressed(self):self.assertTrue(edge_review()['checks']['normalized_stress_not_assumed_bounded_at_outer_edge']);self.assertIn('delta^-3',protocol()['edge'])
    def test_25_exact_half_weight_constants(self):self.assertEqual(edge_review()['half_weight_constants'],['t1^2/2','2'])
    def test_26_no_bounded_normalized_edge_or_endpoint_gap_claim(self):self.assertFalse(FLAGS['normalized_target_or_root_quotient_bounded_on_closed_shell']);self.assertFalse(FLAGS['strict_positive_root_gap_at_zero_stress_endpoints'])
    def test_27_no_core_frame_or_unconditional_promotion(self):self.assertFalse(FLAGS['all_retained_core_background_and_frame_floors_independently_proved']);self.assertFalse(FLAGS['unconditional_actual_pulse_family_or_blowup_verified'])
    def test_28_frozen_record_and_source_bindings(self):
        frozen=json.loads((HERE/'evidence.json').read_text());validate(frozen);compare_records(self.record,frozen)
        self.assertEqual(len(self.record['source_hashes']),14);self.assertGreater(len(self.record['inherited_provenance']),3100)
    def test_29_reject_false_scope_promotions(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_30_reject_critical_formula_mutation(self):
        for k in ('sharp_positive_threshold','edge','Haar_scale'):
            r=deepcopy(self.record);r['protocol'][k]='zero'
            with self.assertRaises(ValueError):validate(r)
    def test_31_reject_missing_history_and_obligations(self):
        for k in ('source_hashes','inherited_provenance','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[k],dict):r[k].pop(next(iter(r[k])))
            else:r[k].pop()
            with self.assertRaises(ValueError):validate(r)
        self.assertTrue(all(parent_review()[1].values()))
    def test_32_reject_failed_boundary_or_gaussian_control(self):
        for key in ('exact_perturbed_cone_and_sharp_obstruction','manufactured_homogeneous_Gaussian_Haar_columns'):
            r=deepcopy(self.record);r['reviews'][key]['checks'][next(iter(r['reviews'][key]['checks']))]=False
            with self.assertRaises(ValueError):validate(r)
    def test_33_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_34_fresh_process_local_import(self):
        code='import sys;sys.path.insert(0,'+repr(str(HERE))+');before=sys.path[:];import audit;assert str(audit.HERE)=='+repr(str(HERE))+';assert before==sys.path;assert "root_bounds" not in sys.modules;assert "covariance_bounds" in sys.modules'
        subprocess.run([sys.executable,'-c',code],check=True,capture_output=True,text=True)
    def test_35_roundtrip_and_all_gates(self):self.assertTrue(all(self.record['gates'].values()));validate(json.loads(json.dumps(self.record)))
    def test_36_negative_numerator_uses_lower_determinant(self):
        gamma,eta=F(1,2),F(2,5)
        bounds=cone_bounds(gamma,eta)
        witness=solve(matrix((0,eta,0,-eta)),(1,1-gamma))[0]
        self.assertLess(witness,bounds['positive_numerator']/bounds['det_upper'])
        self.assertLessEqual(bounds['root_lower'],witness)
        self.assertTrue(corner_review()['checks']['negative_numerator_lower_bound_uses_smaller_determinant'])
        self.assertEqual(corner_review()['negative_numerator_corner_target_cases'],32)
        for e in product((-eta,eta),repeat=4):
            for W in (-(1-gamma),1-gamma):
                self.assertTrue(all(bounds['root_lower']<=q<=bounds['root_upper'] for q in solve(matrix(e),(1,W))))
if __name__=='__main__':unittest.main()
