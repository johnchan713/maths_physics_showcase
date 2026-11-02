import unittest,json
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import sympy as s
import audit,q111
text=(audit.HERE.parent/'core_current_head_q110_radial_tail_request/test_audit.py').read_text()
ns=dict(q110=q111.radial,audit=audit,Path=Path,Q=Q)
exec(text[text.index('def setup('):text.index('\nclass Tests(')],ns);setup=ns['setup']

class Tests(unittest.TestCase):
    def test_forward_majorant_dominates_exact_frozen_recurrence(self):
        frozen=q111.radial.axial.angular.angular
        for H,lo in ((Q(1,5),Q(1,7)),(Q(8),Q(1)),(Q(2),Q(3))):
            tau=min(Q(1),lo)/(16*(1+H));x=tau/lo
            A=4*tau+H*((1-x)**-3-1)
            self.assertLess(A,Q(1,2))
            for n in range(13):
                exact=frozen.native_factor(n,H,lo)
                bound=q111.materialize(q111.forward_recipe(n,H,lo))
                self.assertGreaterEqual(bound,exact)
        for j in range(1,50):self.assertLessEqual(comb(j+2,2),3**j)
    def test_reverse_majorant_dominates_exact_exponential_recurrence(self):
        root=audit.HERE.parent/'core_weighted_reverse_chart/bounds.py'
        import importlib.util
        spec=importlib.util.spec_from_file_location('exact_reverse',root);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        for a,b in ((Q(1,7),Q(3)),(Q(1),Q(8)),(Q(2),Q(4))):
            chart=m.ReverseChart(a,b);B=max(Q(1),b*b,1/a);tau=min(Q(1),a)/(16*(1+B))
            self.assertLess(4*tau+3*B*tau/(1-tau/a),Q(1,2))
            for n in range(13):
                self.assertGreaterEqual(q111.materialize(q111.reverse_recipe(n,a,b)),chart.transfer(n).evaluate(1))
    def test_power_chart_actual_polynomials_both_directions(self):
        U,R,y=s.symbols('U R y',positive=True)
        # R=U^2 corresponds to d=1/2; all other coordinates are unchanged.
        f=R**4*(1+y+y*y);g=f.subs(R,U**2)
        for n in range(6):
            nr=sum(abs(s.diff(f,R,i,y,j).subs({R:4,y:1}))/factorial(i)/factorial(j) for i in range(n+1) for j in range(n+1-i))
            nu=sum(abs(s.diff(g,U,i,y,j).subs({U:2,y:1}))/factorial(i)/factorial(j) for i in range(n+1) for j in range(n+1-i))
            self.assertLessEqual(nu,q111.materialize(q111.forward_recipe(n,4,2))*nr)
            # Independently check U=R^2 in the reverse direction, at R=2.
            raw_u=sum(abs(s.diff(U**4*(1+y+y*y),U,i,y,j).subs({U:4,y:1})) for i in range(n+1) for j in range(n+1-i))
            out_r=sum(abs(s.diff(R**8*(1+y+y*y),R,i,y,j).subs({R:2,y:1}))/factorial(i)/factorial(j) for i in range(n+1) for j in range(n+1-i))
            self.assertLessEqual(out_r,q111.materialize(q111.reverse_recipe(n,1,3))*raw_u)
    def test_same_hull_enclosure_and_large_deferred_plan(self):
        obj,b=setup();bridge=q111.ChartBridge(obj);out=bridge.plan(0)
        self.assertIs(out['source_state'],obj.state)
        self.assertEqual(out['source_chart']['factorial_argument'],100005)
        self.assertEqual(out['additional_raw_source_factorial'],100005)
        self.assertEqual(out['source_factorial_multiplicity'],2)
        self.assertEqual(out['output_chart']['factorial_argument'],0)
        self.assertFalse(out['source_chart']['materialized']);self.assertFalse(out['radial_tail_certified'])
        with self.assertRaises(ValueError):q111.materialize(out['source_chart'])
        for a,bv in ((Q(1,4),Q(3)),(Q(2),Q(7))):
            for d in (Q(1,2),Q(1),Q(2)):
                h=q111.hull(a,bv,d)
                for x in (a,bv):
                    self.assertLessEqual(s.Rational(h['Umin']),s.Rational(x)**s.Rational(d))
                    self.assertGreaterEqual(s.Rational(h['Umax']),s.Rational(x)**s.Rational(d))
        obj.state=object()
        with self.assertRaises(ValueError):bridge.plan(0)
    def test_recipe_retains_old_kernel_gates_and_both_factorials(self):
        obj,b=setup(9);bridge=q111.ChartBridge(obj)
        out=bridge.recipe(0,b,geometry=dict(width=Q(1,4),gap=Q(1,8),L=2),kappa=Q(1,5))
        self.assertEqual(out['raw_tail']['source_factorial'],9)
        self.assertEqual(out['charts']['source_chart']['factorial_argument'],9)
        self.assertFalse(out['raw_tail']['materialized_actual_input_bound'])
        self.assertFalse(out['full_cone_or_actual_pulse_certified'])
        low,lb=setup()
        with self.assertRaises(ValueError):q111.ChartBridge(low).recipe(0,lb,geometry=dict(width=Q(1,4),gap=Q(1,8),L=2))
    def test_invalid_chart_inputs(self):
        for a,b,d in ((0,1,1),(2,1,1),(1,2,Q(1,3)),(1,2,Q(5,2))):
            with self.assertRaises(ValueError):q111.hull(a,b,d)
        for n in (-1,True):
            with self.assertRaises(ValueError):q111.forward_recipe(n,1,1)
        with self.assertRaises(ValueError):q111.forward_recipe(0,1,0)
        with self.assertRaises(TypeError):q111.ChartBridge(None)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
