import unittest
from fractions import Fraction as Q
import json
import ast
from copy import deepcopy
import audit

class ParameterTests(unittest.TestCase):
    def test_selected_parameter_is_outside_constructed_upper_bound(self):
        self.assertLess(Q(1,513),Q(1,200))
    def test_published_sample_reproduces(self):
        self.assertEqual(audit.cross_exponent(Q(1,200),12817),Q(20030773,20000000))
    def test_actual_parameter_upper_bound_fails_stage(self):
        self.assertLess(audit.cross_exponent(Q(1,513),12817),0)
    def test_monotonicity_for_frozen_stage(self):
        slope=(Q(1,2)-3*Q(1,100000))*12818/2-4
        self.assertGreater(slope,0)
    def test_stage_floor_is_sharp(self):
        for h in (Q(1,200),Q(1,1000),Q(1,10**30)):
            for m,N in ((0,0),(2,3),(9,11)):
                j=audit.stage_floor(h,m,N)
                self.assertGreaterEqual(audit.cross_exponent(h,j,m),N+1)
                self.assertLess(audit.cross_exponent(h,j-1,m),N+1)
    def test_formula_matches_frozen_provider_functions(self):
        # Execute only actual pure arithmetic definitions, without loading the
        # expensive unrelated research imports; compare against independent formula.
        path=audit.RESULTS/'core_complete_potential_diagonal_join/diagonal.py'
        tree=ast.parse(path.read_text())
        names={'index','tuple_loss','cut_loss','total_loss'}
        body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
        ns={'Q':Q,'parameters':lambda h:{'A':Q(1,2)+h}}
        exec(compile(ast.Module(body=body,type_ignores=[]),str(path),'exec'),ns)
        for h in (Q(1,200),Q(1,513)):
            for m in (0,2,9):
                loss=ns['cut_loss'](m+2,h)+ns['total_loss'](m+2,h)
                self.assertEqual(loss,15+4*h+4*m)
    def test_invalid_domains(self):
        for h in (0,-1,Q(1,100)):
            with self.assertRaises(ValueError):audit.cross_exponent(h,1)
    def test_frozen_evidence_and_scope(self):
        record=audit.generate()
        self.assertEqual(record,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(record['same_parameter_join_valid'])
        self.assertFalse(record['continuum_theorem_disproved'])
        self.assertFalse(record['endpoint_bound_proved'])

if __name__=='__main__':unittest.main()
