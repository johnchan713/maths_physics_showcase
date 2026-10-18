"""Signed moments, exact boundary controls, source conventions and scope guards."""
from copy import deepcopy
from fractions import Fraction as F
import importlib.util,json,tempfile,unittest
from pathlib import Path
import audit
import moments as m

class SignedPressureMomentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def z(self,R=F(1),angle=F(2,5),nu=F(1)):return m.identities(1,angle,F(1,200),nu,R)
    def test_all_80_exact_cases(self):
        rows=self.r['exact_full_and_truncated_support']['rows'];self.assertEqual(len(rows),80);self.assertTrue(all(all(r['checks'].values()) for r in rows))
    def test_320_scalar_jet_identities(self):self.assertEqual(self.r['exact_full_and_truncated_support']['exact_integrated_scalar_jet_identities'],320)
    def test_pointwise_conservative_flux(self):self.assertTrue(self.z()['pointwise'])
    def test_unweighted_full_support_identity(self):
        z=self.z();self.assertEqual(z['direct'][0].pair(),z['reduced'][0].pair())
    def test_weighted_full_support_identity(self):
        z=self.z();self.assertEqual(z['direct'][1].pair(),z['reduced'][1].pair())
    def test_full_support_boundary_jets_vanish(self):
        z=self.z();self.assertTrue(z['support_endpoint']);self.assertTrue(all(j.pair()==(0,0) for j in z['boundary']))
    def test_truncated_boundary_correction_unweighted(self):
        z=self.z(F(1,3));self.assertEqual(z['direct'][0].pair(),(z['reduced'][0]+z['boundary'][0]).pair());self.assertNotEqual(z['boundary'][0].pair(),(0,0))
    def test_truncated_boundary_correction_weighted(self):
        z=self.z(F(1,3));self.assertEqual(z['direct'][1].pair(),(z['reduced'][1]+z['boundary'][1]).pair());self.assertNotEqual(z['boundary'][1].pair(),(0,0))
    def test_axis_viscosity_not_erased(self):self.assertTrue(self.z()['omission_axis'])
    def test_axis_viscosity_linear_in_nu(self):self.assertEqual(self.z(nu=F(2))['axis'].pair(),(2*self.z()['axis']).pair())
    def test_inviscid_axis_term_zero(self):self.assertEqual(self.z(nu=F(0))['axis'].pair(),(0,0))
    def test_full_weighted_total_independent_of_viscosity(self):self.assertEqual(self.z(nu=F(1))['direct'][1].pair(),self.z(nu=F(3))['direct'][1].pair())
    def test_unweighted_viscosity_total_is_axis_flux(self):
        a=self.z(nu=F(1));b=self.z(nu=F(0));self.assertEqual((a['direct'][0]-b['direct'][0]).pair(),a['axis'].pair())
    def test_premature_future_cutoff_detected(self):self.assertTrue(self.z(F(1,3))['omission_future'])
    def test_both_angular_endpoints(self):
        for e in (F(-1),F(1)):self.assertTrue(self.z(angle=e)['exact'])
    def test_midplane_still_has_axis_flux(self):self.assertTrue(self.z(angle=F(0))['omission_axis'])
    def test_original_source_18_independent_comparisons(self):self.assertEqual(self.r['original_source_comparison']['samples'],18);self.assertTrue(all(self.r['original_source_comparison']['checks'].values()))
    def test_reciprocal_L_first_jet(self):self.assertEqual(m.source_bounds(F(1,200),1,7,3,5)['reciprocal_L_1'],F(10100,9801))
    def test_reciprocal_L_second_jet(self):self.assertTrue(self.r['norm_estimates']['checks']['all_reciprocal_and_average_bounds'])
    def test_uniform_h_interval_not_only_samples(self):self.assertTrue(self.r['norm_estimates']['checks']['uniform_h_interval_follows_from_monotone_positive_majorant'])
    def test_exact_unweighted_bound_fixture(self):self.assertEqual(m.source_bounds(F(1,200),1,7,3,5)['I0'],F(6773825,6534))
    def test_exact_weighted_bound_fixture(self):self.assertEqual(m.source_bounds(F(1,200),1,7,3,5)['I1'],F(145233305,39204))
    def test_signed_bound_axis_viscosity_cost(self):self.assertTrue(self.r['norm_estimates']['checks']['axis_viscosity_cost_is_not_erased'])
    def test_signed_bound_weighted_viscosity_cancels(self):self.assertTrue(self.r['norm_estimates']['checks']['first_weighted_total_has_no_viscosity_boundary'])
    def test_signed_bound_does_not_bound_absolute_integral(self):self.assertTrue(self.r['signed_not_absolute_control']['checks']['absolute_full_integral_is_strictly_larger_than_Q1'])
    def test_auxiliary_viscosity_control_keeps_signed_identity(self):self.assertTrue(self.r['signed_not_absolute_control']['checks']['signed_full_moment_obeys_Q1']);self.assertTrue(self.r['signed_not_absolute_control']['manufactured_auxiliary_viscosity_not_actual_parameters'])
    def test_zero_global_axial_norm_gives_zero_sources(self):
        r=m.from_global_axial_norm(F(1,200),1,7,0);self.assertEqual(r['I0'],0);self.assertEqual(r['I1'],0)
    def test_single_global_axial_norm_controls_radial_average(self):
        r=m.from_global_axial_norm(F(1,200),1,7,3);self.assertEqual(r['u2'],3);self.assertEqual(r['v2'],108)
    def test_angular_loss_factors(self):
        r=self.r['norm_estimates']['checks'];self.assertTrue(r['first_derivative_Cauchy_norm_factor_two']);self.assertTrue(r['second_derivative_norm_factor_three'])
    def test_signed_pressure_source_half_in_m3(self):self.assertEqual(m.repair_moment_bounds(2,3,5,7,11,13)[2],F(131,2))
    def test_weighted_source_half_in_m5(self):self.assertEqual(m.repair_moment_bounds(2,3,5,7,11,13)[4],F(301,2))
    def test_repair_full_memory_not_zeroed(self):
        a=m.repair_moment_bounds(2,3,5,7,11,13);b=m.repair_moment_bounds(2,3,5,7,0,0)
        self.assertEqual(a[2]-b[2],F(11,2));self.assertEqual(a[4]-b[4],F(13,2))
    def test_same_parent_linear_repair_interface(self):
        p=audit.HERE.parent/'first_seminorm_bound/seminorms.py';s=importlib.util.spec_from_file_location('signed_moment_parent_interface',p);w=importlib.util.module_from_spec(s);s.loader.exec_module(w)
        mu=m.repair_moment_bounds(2,3,5,7,11,13);self.assertEqual(mu,w.moment_bounds(2,3,5,7,(11,13)));self.assertGreater(w.repair_bound(mu,2,3,F(1,1000))['coefficient'],0)
    def test_regular_axis_polynomial_quotient(self):
        p=m.P({1:2,2:3});self.assertEqual(p.over_X().value(F(0)).pair(),(2,0))
    def test_singular_axis_quotient_rejected(self):
        with self.assertRaises(ValueError):m.P({0:1}).over_X()
    def test_invalid_source_parameters(self):
        for values in ((0,1,1,1,1),(F(1,10),1,1,1,1),(F(1,200),-1,1,1,1),(F(1,200),1,0,1,1),(F(1,200),1,1,-1,1)):
            with self.assertRaises(ValueError):m.source_bounds(*values)
    def test_invalid_global_axial_norm(self):
        with self.assertRaises(ValueError):m.from_global_axial_norm(F(1,200),1,1,-1)
    def test_invalid_repair_inputs(self):
        for args in ((0,1,1,1,1,1),(1,1,1,1,-1,1)):
            with self.assertRaises(ValueError):m.repair_moment_bounds(*args)
    def test_exact_inputs_refuse_float_and_bool(self):
        for x in (True,.005):
            with self.assertRaises(TypeError):m.source_bounds(x,1,1,1,1)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'a70349e4b48b41c08033921824951ce286766195');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_bytes_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_no_scope_promotion(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_global_angular_obligation_not_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_integer_boolean_promotion_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_axis_term_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['exact_full_and_truncated_support']['rows'][0]['axis']=['0','0']
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
