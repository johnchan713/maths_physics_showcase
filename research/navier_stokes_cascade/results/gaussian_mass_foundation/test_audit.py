"""Independent Gaussian/Haar controls, negative cases and frozen-scope guards."""
from copy import deepcopy
from fractions import Fraction as F
import json
from math import exp,isfinite
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,main
from gaussian_bounds import curvature,root_interval,rate_curvature_interval,barrier,logs,square_budget,curvature_review,barrier_review,tail_review,degeneracy_review
from gaussian_checks import rate,log_envelope,cutoff,integral_review,ode_mass_review,exact_haar_review


class GaussianMassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()
    def test_current_parent_and_protocol(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
        self.assertEqual(protocol()['parent_commit'],'a70298fb0db13af479b2547b52e5bf1454b99e6d')
    def test_exact_rational_curvature(self):
        c,C=curvature(1,F(1,2),2);self.assertEqual(c,F(1,4)*(F(16,169)+F(1,2)))
        self.assertEqual(C,F(27,5));self.assertGreater(c,0)
    def test_zero_tilt_not_accepted(self):
        for u in (0,-1):
            with self.assertRaises(ValueError):curvature(u,1,2)
    def test_degenerate_spectral_bounds_rejected(self):
        for lo,hi in ((0,1),(2,1),(-1,2)):
            with self.assertRaises(ValueError):curvature(1,lo,hi)
    def test_exact_square_enclosure(self):self.assertEqual(root_interval(F(9,16)),(F(3,4),F(3,4)))
    def test_irrational_square_enclosure(self):
        lo,hi=root_interval(2);self.assertLess(lo*lo,2);self.assertGreater(hi*hi,2)
        self.assertEqual(hi-lo,F(1,2**80))
    def test_invalid_square_enclosure(self):
        for x,bits in ((-1,80),(1,0),(1,F(3,2))):
            with self.assertRaises(ValueError):root_interval(x,bits)
    def test_whole_pulse_curvature_control_cases(self):
        r=curvature_review();self.assertEqual(r['rational_enclosure_cases'],75);self.assertTrue(all(r['checks'].values()))
    def test_interval_endpoints_are_included(self):
        for a in (F(1,2),F(3,2)):
            l,h=rate_curvature_interval(1,a,1);self.assertGreater(l,0);self.assertGreaterEqual(h,l)
    def test_outside_pulse_rate_rejected(self):
        with self.assertRaises(ValueError):rate_curvature_interval(1,2,1)
    def test_center_rate_and_log_normalization(self):
        self.assertAlmostEqual(rate(50,100,1,1),0,places=13);self.assertEqual(log_envelope(50,100,1,1),0)
    def test_sign_of_rate_either_side_of_center(self):
        self.assertGreater(rate(25,100,1,1),0);self.assertLess(rate(75,100,1,1),0)
    def test_independent_closed_integral_and_derivative(self):self.assertTrue(all(integral_review()['checks'].values()))
    def test_Riccati_barrier_sharp_padded_slack(self):
        r=barrier(1,1,F(1,4),8);self.assertEqual(r['rho'],F(1,8));self.assertEqual(r['outward'],r['negative_bound'])
    def test_unaccepted_Riccati_length_rejected(self):
        with self.assertRaises(ValueError):barrier(1,1,0,7)
    def test_negative_error_bound_rejected(self):
        with self.assertRaises(ValueError):barrier(1,1,-1,8)
    def test_total_logarithmic_amplitude_budget(self):self.assertEqual(barrier(1,1,F(1,4),16)['D'],F(11,8))
    def test_manufactured_noncommuting_frame_control(self):self.assertTrue(ode_mass_review()['checks']['noncommuting_frame_fixture'])
    def test_actual_cutoff_Haar_integral(self):
        r=ode_mass_review();self.assertEqual(len(r['rows']),8)
        self.assertTrue(r['checks']['actual_cutoff_Haar_integrals_between_explicit_bounds'])
        self.assertTrue(r['checks']['actual_first_absolute_moment_concentration'])
    def test_original_cutoff_plateau_and_support(self):
        self.assertEqual(cutoff(.5),1);self.assertEqual(cutoff(.3),1)
        self.assertEqual(cutoff(0),0);self.assertEqual(cutoff(1),0)
        self.assertGreater(cutoff(.25),0);self.assertLess(cutoff(.25),1)
    def test_Haar_radius_scaling_and_angular_half(self):self.assertTrue(all(exact_haar_review()['checks'].values()))
    def test_log_mass_floor_survives_underflow(self):
        r=logs(1,F(1,2),2,10000,10000);self.assertTrue(isfinite(r['log_jminus']));self.assertEqual(exp(r['log_jminus']),0)
    def test_tilt_limit_has_no_uniform_Gaussian_curvature(self):self.assertTrue(all(degeneracy_review()['checks'].values()))
    def test_exact_tail_rate_keeps_frozen_radius(self):
        r=tail_review();c,_=curvature(1,F(1,2),2);self.assertEqual(F(r['a0']),c*F(1,6057072)/25)
    def test_square_completion_no_new_first_band(self):
        d,K=square_budget(F(1,10000),17,10,30)
        for ell in (1,10,100,10000):self.assertLessEqual(-F(1,10000)*ell*ell+d*ell,K)
    def test_negative_physical_loss_budget(self):self.assertEqual(square_budget(F(1,3),0,-5,1),(0,0))
    def test_invalid_tail_budget_rejected(self):
        for args in ((0,1,1,1),(1,-1,1,1),(1,1,1,-1)):
            with self.assertRaises(ValueError):square_budget(*args)
    def test_pressure_memory_and_tail_source_boundary_kept(self):
        self.assertIn('Pi1/Pi2',protocol()['memory']);self.assertIn('never forwarded',protocol()['memory'])
        self.assertIn('C^-20',protocol()['fixed_choices']);self.assertIn('1/3',protocol()['support'])
    def test_frozen_evidence_and_all_bindings(self):
        frozen=json.loads((HERE/'evidence.json').read_text());validate(frozen);compare_records(self.record,frozen)
        self.assertEqual(len(frozen['source_hashes']),18);self.assertGreaterEqual(frozen['inherited_provenance']['records'],76)
    def test_scope_flags_cannot_be_promoted(self):
        for key,value in FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):validate(r)
    def test_failed_numerical_control_rejected(self):
        r=deepcopy(self.record);r['reviews']['noncommuting_frame_and_actual_cutoff_mass']['checks']['actual_cutoff_Haar_integrals_between_explicit_bounds']=False
        with self.assertRaises(ValueError):validate(r)
    def test_changed_curvature_and_dropped_limit_rejected(self):
        for key in ('protocol','remaining_obligations'):
            r=deepcopy(self.record)
            if key=='protocol':r[key]['curvature']='uniform at zero tilt'
            else:r[key]=[]
            with self.assertRaises(ValueError):validate(r)
    def test_history_binding_deletion_rejected(self):
        r=deepcopy(self.record);r['inherited_provenance']['bindings'].pop(next(iter(r['inherited_provenance']['bindings'])))
        with self.assertRaises(ValueError):validate(r)
    def test_float_tolerance_does_not_allow_changed_algebra(self):
        r=deepcopy(self.record);r['reviews']['whole_pulse_curvature']['rational_enclosure_cases']=74
        with self.assertRaises(ValueError):validate(r)
    def test_frozen_output_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_json_roundtrip(self):validate(json.loads(json.dumps(self.record)))
if __name__=='__main__':unittest.main()
