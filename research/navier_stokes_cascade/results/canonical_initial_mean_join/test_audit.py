"""Canonical complete real covariance, pressure and initial-only scope checks."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,main
from real_covariance import jet,dot,outer,add,scale,sub,covariance,phase_average,complete_review,variance_case,variance_review
from mean_checks import pressure_data,pressure_review,order_review,nonflat_review,recompute_review


class CanonicalInitialMeanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()
    def test_current_parent_and_frozen_protocol(self):
        self.assertEqual(protocol()['parent_commit'],'d3c4dd1b00337d54d446888355802a9d8a824bf8')
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
    def test_independent_exact_phase_quadrature(self):
        r=complete_review();self.assertEqual(len(r['rows']),18);self.assertTrue(all(r['checks'].values()))
    def test_old_and_new_amplitudes_share_transverse_plane(self):
        for eta in (F(0),F(1,7)):
            n,A,B,fp=jet(F(1,5),eta);self.assertEqual(dot(n,A),0);self.assertEqual(dot(n,B),fp)
    def test_full_curl_square_is_not_principal_only(self):
        _,A,B,_=jet(F(1,5));Q=covariance(A,B,11)
        self.assertNotEqual(Q,scale(outer(A,A),F(1,2)))
        self.assertEqual(Q,phase_average(A,B,11))
    def test_curl_square_k_squared_scaling(self):
        _,A,B,_=jet(F(1,5));principal=scale(outer(A,A),F(1,2))
        self.assertEqual(scale(sub(covariance(A,B,3),principal),9),scale(sub(covariance(A,B,11),principal),121))
    def test_all_six_covariance_changes_retained(self):self.assertTrue(complete_review()['checks']['all_six_tensor_change_identity'])
    def test_generic_complex_sine_cross_not_cancelled(self):
        self.assertTrue(complete_review()['checks']['generic_complex_phase_shift_has_nonzero_linear_sine_cross'])
        self.assertFalse(FLAGS['real_primary_cancellation_applies_to_generic_complex_later_wave'])
    def test_positive_rebuilt_variance_weights(self):
        r=variance_case(F(1,20));self.assertGreater(min(r['weights']),0)
        self.assertEqual(sum(r['weights']),1)
    def test_degenerate_nonpositive_weight_boundary_rejected(self):
        for eta in (F(-1,2),F(-1)):
            with self.assertRaises(ValueError):variance_case(eta)
    def test_same_two_fluxes_leave_nonzero_variance(self):
        r=variance_case(F(1,20));self.assertEqual(r['Wa'][0][1:3],r['Wr'][0][1:3])
        self.assertEqual(r['delta'],r['expected']);self.assertNotEqual(r['delta'][2][2],0)
    def test_zero_variance_change_at_zero_eta(self):
        self.assertEqual(variance_case(0)['delta'],[[0]*3 for _ in range(3)])
    def test_stale_roots_do_not_even_match_two_fluxes(self):self.assertTrue(variance_review()['checks']['stale_weights_also_break_selected_flux'])
    def test_full_pressure_and_moment_identities(self):self.assertTrue(all(pressure_review()['checks'].values()))
    def test_pressure_compact_and_bump_unit_mass(self):
        d=pressure_data(F(1,20));self.assertEqual(d['rho'].integral(),1)
        self.assertEqual(d['pm'].at(0,1),0);self.assertEqual(d['pm'].at(0,2),0)
    def test_axial_pressure_moment_not_zero(self):
        d=pressure_data(F(1,20));self.assertTrue(d['moment']);self.assertTrue(d['moment'].diff(1))
    def test_omitted_moving_bump_derivative_is_nonzero(self):
        d=pressure_data(F(1,20));self.assertTrue(d['crho'].diff(1)*d['P'])
    def test_all_three_fresh_defects_nonzero(self):
        d=pressure_data(F(1,20));self.assertTrue(all(d[n] for n in ('P','Jtheta','Jz')))
    def test_fresh_five_rows_cancel_current_not_old_targets(self):
        r=pressure_review();self.assertTrue(r['checks']['fixed_mean_matrix_cancels_recomputed_linear_targets'])
        self.assertTrue(r['checks']['same_two_principal_fluxes_can_leave_three_nonzero_new_defects'])
    def test_real_curl_error_order_is_quadratic(self):
        r=order_review();self.assertEqual(F(r['orders']['real_curl_covariance_square']),2-2*F(r['kappa']))
    def test_initial_mean_viscosity_loss_kept(self):
        r=order_review();self.assertEqual(F(r['orders']['final_tangential_mean']),2-3*F(r['kappa']))
        self.assertTrue(all(r['checks'].values()))
    def test_initial_defects_S2_not_zero(self):
        r=order_review();self.assertEqual(r['orders']['final_defects'],'2')
        self.assertFalse(FLAGS['nonlinear_defects_identically_zero_after_five_linear_rows'])
    def test_exact_polynomial_not_epsilon_gain_controls(self):
        r=nonflat_review();self.assertEqual(len(r['rows']),6);self.assertTrue(all(r['checks'].values()))
        self.assertTrue(all(x['log2_ratio_lower_bound']>0 for x in r['rows']))
    def test_initial_Haar_and_nonlinear_recompute(self):self.assertTrue(all(recompute_review()['checks'].values()))
    def test_later_nonzero_Haar_state_does_not_get_initial_cancellation(self):
        self.assertTrue(recompute_review()['checks']['nonzero_old_mean_breaks_sharp_reduction'])
    def test_protocol_keeps_operation_order(self):
        self.assertIn('recompute',protocol()['order']);self.assertIn('current targets',protocol()['order'])
        self.assertIn('canonical',protocol()['construction'])
    def test_full_pressure_and_axial_velocity_remain(self):
        self.assertIn('Pi1/Pi2',protocol()['memory']);self.assertIn('axial velocity',protocol()['memory'])
        self.assertIn('whole partial_Z',protocol()['memory'])
    def test_fixed_width_actual_flat_fields_and_hull_remain(self):
        self.assertIn('C^-20',protocol()['construction']);self.assertIn('A0/A1',protocol()['support'])
        self.assertIn('I4 potential hull',protocol()['support'])
    def test_frozen_evidence_and_sources(self):
        old=json.loads((HERE/'evidence.json').read_text());validate(old);compare_records(self.record,old)
        self.assertEqual(len(old['source_hashes']),18);self.assertGreaterEqual(old['inherited_provenance']['records'],79)
    def test_scope_promotions_rejected(self):
        for key,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):validate(r)
    def test_erased_history_and_limits_rejected(self):
        for key in ('inherited_provenance','remaining_obligations'):
            r=deepcopy(self.record)
            if key=='inherited_provenance':r[key]['bindings'].pop(next(iter(r[key]['bindings'])))
            else:r[key]=[]
            with self.assertRaises(ValueError):validate(r)
    def test_exact_record_corruption_rejected(self):
        r=deepcopy(self.record);r['reviews']['fresh_pressure_moments_and_targets']['rows'][0]['fresh_defects'][0]='0'
        with self.assertRaises(ValueError):validate(r)
    def test_protocol_and_failed_gate_rejected(self):
        r=deepcopy(self.record);r['protocol']['real_covariance']='principal only'
        with self.assertRaises(ValueError):validate(r)
        r=deepcopy(self.record);r['gates']['canonical_complete_real_curl']=False
        with self.assertRaises(ValueError):validate(r)
    def test_frozen_output_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_json_roundtrip(self):validate(json.loads(json.dumps(self.record)))
if __name__=='__main__':unittest.main()
