"""Exact row, pressure, exterior-memory, provenance and scope regressions."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from identities import annular_rows,annular_inverse,compact_rows,compact_inverse,reviews
from audit import HERE,FLAGS,REMAINING,compare_records,dependency_review,generate,main,protocol,validate


class GlobalPressureMomentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=generate()
    def test_01_parent(self):self.assertEqual(protocol()['parent_commit'],'29b1f5f281650f5fae082d5d4a90a1f0d66e9ae3')
    def test_02_protocol(self):self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
    def test_03_annular_transform_inverse(self):
        x=tuple(map(F,('2','3','5','7','11')));self.assertEqual(annular_inverse(annular_rows(x,F(2,5)),F(2,5)),x)
    def test_04_annular_zero_equivalence(self):self.assertEqual(annular_inverse((F(0),)*5,F(-1)),(F(0),)*5)
    def test_05_compact_transform_inverse(self):
        x=tuple(map(F,('2','3','5','7','11')));self.assertEqual(compact_inverse(compact_rows(x,F(1,10**30)),F(1,10**30)),x)
    def test_06_zero_lambda_rejected(self):
        with self.assertRaises(ValueError):compact_rows((F(0),)*5,F(0))
    def test_07_all_exact_identity_reviews(self):self.assertTrue(all(all(v['checks'].values()) for v in reviews().values()))
    def test_08_pressure_negative_control(self):self.assertTrue(reviews()['same_gauge_complete_pressure_ledger']['checks']['omitted_pressure_repair_breaks_backward_identity'])
    def test_09_pressure_eta_derivative_negative_control(self):self.assertTrue(reviews()['zero_defect_exterior_stress_join']['checks']['pressure_eta_derivative_retained'])
    def test_10_dependencies_reproduced(self):
        for n in ('matching_root_foundation','global_moment_stress','regional_cone_join'):
            r,h=dependency_review(n);self.assertTrue(all(r['gates'].values()) and all(h.values()))
    def test_11_current_scope_closes_only_current_pressure_moment_obligation(self):
        self.assertTrue(FLAGS['complete_core_backward_pressure_same_gauge_trace_discharged'])
        self.assertFalse(FLAGS['unconditional_uniform_actual_pulse_estimates_verified'])
        self.assertFalse(FLAGS['unconditional_complete_residual_or_blowup_verified'])
        self.assertEqual(len(REMAINING),4)
    def test_12_frozen_record_and_bindings(self):
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.r,f);validate(json.loads(json.dumps(self.r)))
        self.assertEqual(len(self.r['source_hashes']),12);self.assertGreater(len(self.r['inherited_provenance']),1000)
    def test_13_scope_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.r);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_14_mutations_rejected(self):
        for key in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.r)
            if isinstance(r[key],dict):r[key].pop(next(iter(r[key])))
            else:r[key].pop()
            with self.assertRaises(ValueError):validate(r)
        r=deepcopy(self.r);r['reviews']['same_gauge_complete_pressure_ledger']['checks']['same_gauge_gives_zero_total_pressure_defect']=False
        with self.assertRaises(ValueError):validate(r)
    def test_15_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
if __name__=='__main__':unittest.main()
