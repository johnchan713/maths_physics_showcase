"""Canonical entrance, full physical equation and root stability regressions."""
from copy import deepcopy
from fractions import Fraction as F
import json
from math import isfinite
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,main
from canonical_bounds import parent_budget,canonical,comparison_logs,budget_review,matrix_review,acceptance_review,mul,sub,inv,norminf,logadd
from canonical_checks import physical_review,cutoff,logP,coefficients,phase_module


class CanonicalPrimaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()
    def test_parent_and_protocol(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
        self.assertEqual(protocol()['parent_commit'],'817edc46cb8ba64db62f3e7c0287c45522ef5eab')
    def test_reference_not_actual_frame_error(self):
        p=parent_budget().inputs();b=parent_budget().budgets(p);r=canonical(p)
        self.assertEqual(r['CE'],b['CEref']);self.assertNotEqual(r['CE'],b['CE'])
    def test_canonical_column_formula(self):self.assertTrue(all(budget_review()['checks'].values()))
    def test_pressure_difference_keeps_denominators(self):
        p=parent_budget().inputs();r=canonical(p)
        smaller=canonical(parent_budget().inputs(bmin=p['bmin']/2))
        self.assertGreater(smaller['pressure_difference_multiplier'],r['pressure_difference_multiplier'])
    def test_Riccati_acceptance_and_outward_sign(self):self.assertTrue(all(acceptance_review()['checks'].values()))
    def test_zero_margin_rejected(self):
        for margin in (0,-1,2):
            with self.assertRaises(ValueError):comparison_logs(parent_budget().inputs(),margin)
    def test_invalid_Haar_floors_rejected(self):
        for lo,hi in ((0,1),(2,1)):
            with self.assertRaises(ValueError):comparison_logs(parent_budget().inputs(),Jmin=lo,Jmax=hi)
    def test_logarithmic_large_constants_stay_finite(self):
        r=comparison_logs(parent_budget().inputs());self.assertTrue(all(isfinite(x) for x in r.values()))
        self.assertGreater(r['log_relative_root_acceptance_length'],1000)
    def test_logadd_avoids_overflow(self):self.assertAlmostEqual(logadd(1000000.,0.),1000000.)
    def test_exact_matrix_cases_and_relative_stability(self):
        r=matrix_review();self.assertEqual(len(r['rows']),24);self.assertTrue(all(r['checks'].values()))
    def test_relative_comparison_retains_radius(self):
        r=matrix_review();self.assertEqual(len({x['radius'] for x in r['rows']}),2)
        self.assertTrue(r['checks']['relative_comparison_retains_and_cancels_common_physical_Haar_scale'])
    def test_stale_roots_not_identified(self):self.assertTrue(all(x['stale_covariance_defect_nonzero'] for x in matrix_review()['rows']))
    def test_new_principal_roots_recover_target(self):self.assertTrue(matrix_review()['checks']['rebuilding_roots_recovers_exact_principal_covariance'])
    def test_inverse_identity_and_singularity(self):
        a=[[F(2),F(1)],[F(1),F(3)]];self.assertEqual(mul(a,inv(a)),[[1,0],[0,1]])
        with self.assertRaises(ValueError):inv([[F(1),F(1)],[F(1),F(1)]])
    def test_matrix_infinity_norm(self):self.assertEqual(norminf([[F(1),F(-2)],[F(2),F(0)]]),3)
    def test_square_root_relative_bound_via_rational_square(self):
        for R in (F(0),F(1,10),F(1,2)):
            self.assertLessEqual((1-R)**2,1-R)
            self.assertGreaterEqual((1+R)**2,1+R)
    def test_reciprocal_root_safe_bound_via_squares(self):
        for R in (F(0),F(1,10),F(1,2)):
            self.assertGreaterEqual((1+2*R)**2*(1-R),1)
    def test_independent_physical_and_frame_integration(self):
        r=physical_review();self.assertEqual(len(r['rows']),4)
        self.assertTrue(r['checks']['independent_reference_frame_and_physical_ODE_agree'])
    def test_complete_pressure_balance_and_omission_control(self):
        r=physical_review();self.assertTrue(r['checks']['full_moving_normal_pressure_satisfies_principal_balance'])
        self.assertTrue(r['checks']['omitting_nprime_pressure_is_detected'])
    def test_actual_cutoff_Haar_columns_and_rebuilt_roots(self):self.assertTrue(all(physical_review()['checks'].values()))
    def test_diagnostic_root_change_is_not_zero(self):
        self.assertTrue(all(x['relative_positive_root_difference']>0 for x in physical_review()['rows']))
    def test_original_cutoff_plateau_support(self):
        self.assertEqual(cutoff(.5),1);self.assertAlmostEqual(cutoff(.2),cutoff(.8))
        self.assertEqual(cutoff(0),0);self.assertEqual(cutoff(1),0)
    def test_envelope_midpoint_normalization(self):self.assertEqual(logP(20,40,.75,.8),0)
    def test_distinct_canonical_and_inherited_initial_columns(self):
        d=phase_module().phase_data(40*3028536,1,0);c=coefficients(d,0)
        self.assertGreater(float(abs(c['Br'][1,0]-c['Ba'][1,0])),0)
    def test_keeps_pressure_and_fixed_geometry(self):
        self.assertIn('Pi1/Pi2',protocol()['memory']);self.assertIn('C^-20',protocol()['fixed'])
        self.assertIn('angular half',protocol()['Haar'])
    def test_no_complete_covariance_or_map_promotion(self):
        self.assertFalse(FLAGS['complete_curl_covariance_equals_principal_covariance'])
        self.assertFalse(FLAGS['later_generated_mean_and_repair_maps_automatically_identified'])
    def test_frozen_evidence_and_sources(self):
        old=json.loads((HERE/'evidence.json').read_text());validate(old);compare_records(self.record,old)
        self.assertEqual(len(old['source_hashes']),18)
    def test_scope_promotions_rejected(self):
        for key,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):validate(r)
    def test_erased_history_rejected(self):
        r=deepcopy(self.record);r['inherited_provenance']['bindings'].pop(next(iter(r['inherited_provenance']['bindings'])))
        with self.assertRaises(ValueError):validate(r)
    def test_erased_limits_rejected(self):
        r=deepcopy(self.record);r['remaining_obligations']=[]
        with self.assertRaises(ValueError):validate(r)
    def test_protocol_and_exact_matrix_corruption_rejected(self):
        r=deepcopy(self.record);r['reviews']['exact_principal_matrix_root_stability']['rows'][0]['theta']='0'
        with self.assertRaises(ValueError):validate(r)
        r=deepcopy(self.record);r['protocol']['canonical']='same pulse'
        with self.assertRaises(ValueError):validate(r)
    def test_ODE_gate_failure_rejected(self):
        r=deepcopy(self.record);r['reviews']['independent_physical_frame_pressure_and_Haar_controls']['checks']['omitting_nprime_pressure_is_detected']=False
        with self.assertRaises(ValueError):validate(r)
    def test_balance_corruption_not_hidden_by_tolerance(self):
        r=deepcopy(self.record);r['reviews']['independent_physical_frame_pressure_and_Haar_controls']['rows'][0]['complete_pressure_balance_error']=1e-5
        with self.assertRaises(ValueError):validate(r)
    def test_frozen_output_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_json_roundtrip(self):validate(json.loads(json.dumps(self.record)))
if __name__=='__main__':unittest.main()
