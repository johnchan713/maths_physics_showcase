"""Exact algebra, actual map, omitted-term, support and immutable-evidence tests."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,parent_review,main
from support_bounds import reciprocal_amplitude_cost,target_bound,map_review,jet_review,support_review,support_edges
from conservative_checks import J,P,X,z2,assemble,conservative_review,plateau_review
from linear_checks import numerical_review

class CoefficientSupportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=generate()
    def test_01_protocol(self):self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
    def test_02_current_parent(self):self.assertEqual(protocol()['parent_commit'],'135c9da9a41f070ce0d583ba2a70d9a8fe21d490')
    def test_03_linear_continuum_inverse(self):self.assertLess(F(map_review()['E_inverse_upper']),800);self.assertLess(F(map_review()['U_inverse_upper']),2)
    def test_04_reciprocal_value_norm(self):self.assertEqual(reciprocal_amplitude_cost(0,F(7)),F(2,7))
    def test_05_reciprocal_first_norm(self):self.assertEqual(reciprocal_amplitude_cost(1,F(7)),F(4,7))
    def test_06_reciprocal_all_higher_norms(self):
        for k in (2,3,10,1000000):self.assertEqual(reciprocal_amplitude_cost(k,F(7)),F(5,7))
    def test_07_target_retains_both_axial_discrepancies(self):
        z,r=target_bound(8,F(7),F(2),F(3,5),F(1,1000),tuple(map(F,(2,3,5,7,11))))
        self.assertEqual(z[1],4000);self.assertEqual(r,3200000)
    def test_08_invalid_norm_inputs(self):
        for k in (-1,True,1.5):
            with self.assertRaises(ValueError):reciprocal_amplitude_cost(k,F(7))
        with self.assertRaises(ValueError):reciprocal_amplitude_cost(2,F(0))
        with self.assertRaises(ValueError):target_bound(2,F(7),F(2),F(3,5),F(0),(1,)*5)
    def test_09_exact_angular_inverse(self):
        a=J((F(2),F(3),F(5)));self.assertEqual((a*a.inverse()).c,(1,0,0))
    def test_10_exact_radial_primitive(self):
        p=P({2:3,5:7});self.assertEqual(p.primitive().dx().total().c,p.total().c)
    def test_11_general_conservative_totals(self):self.assertEqual(conservative_review()['exact_integrated_identity_count'],48);self.assertTrue(all(conservative_review()['checks'].values()))
    def test_12_previous_viscosity_not_deleted(self):self.assertTrue(any(r['omitted_lower_viscosity_nonzero'] for r in conservative_review()['rows']))
    def test_13_additional_order_seven_exact(self):
        a,b,v=assemble(7,F(1,200),F(1,3),F(2));self.assertEqual(a,b);self.assertTrue(any(v))
    def test_14_shifted_radial_pressure_plateau(self):self.assertEqual(plateau_review()['manufactured_second_pressure_plateau'],'153/10');self.assertTrue(all(plateau_review()['checks'].values()))
    def test_15_first_support_edges(self):self.assertEqual(support_edges(1),dict(velocity='Xvel',pressure='Xv',stress='Xb'))
    def test_16_second_support_edges(self):self.assertEqual(support_edges(2),dict(velocity='Xvel',pressure='Xv',stress='Xv'))
    def test_17_later_common_support(self):
        for n in (3,7,1000000):self.assertEqual(set(support_edges(n).values()),{'Xvel'})
    def test_18_support_geometry_and_pressure(self):self.assertTrue(all(support_review()['checks'].values()));self.assertEqual(support_review()['mean_pressure'][1],'K2(eta)')
    def test_19_actual_linear_integrands(self):self.assertEqual(len(numerical_review()['rows']),8);self.assertTrue(numerical_review()['checks']['independent_ordinary_integrands_close_all_five_moments'])
    def test_20_nonlinear_repair_negative_control(self):self.assertTrue(numerical_review()['checks']['nonlinear_leading_repair_rejected_for_positive_order'])
    def test_21_actual_second_jets(self):self.assertTrue(numerical_review()['checks']['same_linear_inverse_controls_second_eta_jets'])
    def test_22_amplitude_jet_negative_control(self):self.assertTrue(numerical_review()['checks']['freezing_actual_amplitude_derivatives_rejected'])
    def test_23_large_targets_allowed(self):self.assertTrue(numerical_review()['checks']['no_small_target_condition'])
    def test_24_scope_not_promoted(self):
        self.assertTrue(FLAGS['every_finite_coefficient_support_and_reserved_mean_velocity_joined'])
        self.assertFalse(FLAGS['reserved_mean_pressure_equal_to_leading_pressure']);self.assertFalse(FLAGS['unconditional_actual_pulse_bounds_verified'])
    def test_25_frozen_record_and_sources(self):
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.r,f)
        self.assertEqual(len(self.r['source_hashes']),16);self.assertGreater(len(self.r['inherited_provenance']),2171)
    def test_26_reject_mutated_scope(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.r);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_27_reject_erased_obligation_gate_and_source(self):
        for k in ('remaining_obligations','gates','source_hashes','protocol'):
            r=deepcopy(self.r)
            if isinstance(r[k],dict):r[k].pop(next(iter(r[k])))
            else:r[k].pop()
            with self.assertRaises(ValueError):validate(r)
    def test_28_reject_history_and_review_mutation(self):
        r=deepcopy(self.r);r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])))
        with self.assertRaises(ValueError):validate(r)
        r=deepcopy(self.r);r['reviews']['actual_linear_map']['inverse_bound']=1000
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(all(parent_review()[1].values()))
    def test_29_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_30_all_reviews(self):self.assertTrue(all(self.r['gates'].values()));self.assertTrue(all(jet_review()['checks'].values()))
    def test_31_standalone_diagnostic_restores_import_search_path(self):
        code='import sys; sys.path.insert(0,'+repr(str(HERE))+'); before=sys.path[:]; import linear_checks; assert sys.path==before; import audit; assert str(audit.HERE)=='+repr(str(HERE))
        subprocess.run([sys.executable,'-c',code],check=True,capture_output=True,text=True)
if __name__=='__main__':unittest.main()
