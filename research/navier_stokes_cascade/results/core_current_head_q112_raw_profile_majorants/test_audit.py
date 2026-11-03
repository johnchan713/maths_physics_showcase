import unittest,json
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import sympy as s
import audit,q112
text=(audit.HERE.parent/'core_current_head_q110_radial_tail_request/test_audit.py').read_text()
ns=dict(q110=q112.charts.radial,audit=audit,Path=Path,Q=Q)
exec(text[text.index('def setup('):text.index('\nclass Tests(')],ns);setup=ns['setup']

def load(name,path):
    import importlib.util
    spec=importlib.util.spec_from_file_location(name,audit.HERE.parent/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

class Tests(unittest.TestCase):
    def test_power_bounds_dominate_original_raw_falling_factorials(self):
        original=load('power_profile','core_current_head_q44_radial_multipliers/radial.py')
        for a,b in ((Q(1,4),Q(3)),(Q(1),Q(8)),(Q(1,2),Q(3,4))):
            for d in (Q(2,5),Q(1),Q(12,5)):
                g=original.RadialGeometry(d,a,b)
                for n in range(11):
                    exact=g.multipliers(n,e=0)
                    self.assertGreaterEqual(q112.materialize(q112.power_recipe(n,a,b,'f0')),exact['f'])
                    self.assertGreaterEqual(q112.materialize(q112.power_recipe(n,a,b,'k0')),exact['k'])
                    self.assertEqual(exact['g'],1)
    def test_actual_fractional_power_derivatives(self):
        U=s.symbols('U',positive=True)
        for d in (Q(1,2),Q(1),Q(2)):
            f=U**(s.Rational(1,d)-1)/s.Rational(d);k=s.Rational(d)*U**(1-s.Rational(1,d))
            for n in range(7):
                for v in (s.Rational(1,4),s.Rational(1),s.Rational(4)):
                    for field,role in ((f,'f0'),(k,'k0')):
                        raw=abs(s.diff(field,U,n).subs(U,v))
                        self.assertLessEqual(raw,q112.materialize(q112.power_recipe(n,Q(1,4),4,role)))
    def test_closed_cutoff_dominates_full_moving_composition(self):
        original=load('moving_profile','core_current_head_q46_moving_cutoff/moving.py')
        for d,a,b,Xb in ((Q(2,5),Q(1,4),Q(3),Q(8)),(Q(1),Q(1),Q(8),Q(8)),(Q(12,5),Q(1,2),Q(1),Q(2))):
            shape=[Q(j+1,2) for j in range(9)]
            moving=original.Moving(d,a,b,Xb,Q(1,4),shape)
            for n in range(9):
                bound=q112.materialize(q112.cutoff_recipe(n,a,Xb),shape_sum=sum(shape[:n+1]))
                self.assertGreaterEqual(bound,max(moving.raw(n)))
                if n:
                    static=original.Moving(d,a,b,Xb,Q(1,4),shape)
                    static.coords.log_s=lambda j:Q(0)
                    self.assertGreater(max(moving.raw(n)),max(static.raw(n)))
    def test_analytic_majorant_radius_and_product_rule(self):
        for a,Xb in ((Q(1,4),Q(8)),(Q(1),Q(1,2)),(Q(2),Q(10))):
            u=min(Q(1),a**3);t=u/(512*(1+Xb))
            L=6*(t/u)/(1-t/u)+Q(2,3)*64*t/(1-64*t)
            self.assertLess(L,1/(8*(1+Xb)))
            self.assertLess(2*Xb*L,Q(1,4))
        U,Z=s.symbols('U Z',positive=True);d=s.Rational(1,2)
        k=d*U**(1-1/d);chi=U**2*(1+Z+Z*Z)
        for n in range(4):
            B=q112.materialize(q112.power_recipe(n,1,4,'k0'))
            C=max(abs(s.diff(chi,U,i,Z,j).subs({U:2,Z:1})) for i in range(n+2) for j in range(n+2-i))
            actual=max(abs(s.diff(k*s.diff(chi,U),U,i,Z,j).subs({U:2,Z:1})) for i in range(n+1) for j in range(n+1-i))
            self.assertLessEqual(actual,2**n*B*C)
    def test_current_binding_and_deferred_three_factorials(self):
        obj,b=setup();p=q112.KernelProfiles(q112.charts.ChartBridge(obj));out=p.plan(0)
        self.assertIs(out['source_state'],obj.state)
        self.assertEqual(out['source_f0']['factorial_argument'],100005)
        self.assertEqual(out['cutoff']['shape_sum_order'],1)
        self.assertEqual(out['source_factorial_multiplicity_including_f0'],3)
        self.assertFalse(out['actual_A0_gain_certified'])
        with self.assertRaises(ValueError):q112.materialize(out['source_f0'])
        with self.assertRaises(ValueError):q112.materialize(out['cutoff'])
        hi,hb=setup(9);proof=q112.KernelProfiles(q112.charts.ChartBridge(hi))
        r=proof.recipe(0,hb,geometry=dict(width=Q(1,4),gap=Q(1,8),L=2),kappa=Q(1,5))
        self.assertEqual(r['profile_majorants']['source_f0']['factorial_argument'],9)
        self.assertFalse(r['full_cone_or_residual_certified'])
        with self.assertRaises(ValueError):p.recipe(0,b,geometry=dict(width=Q(1,4),gap=Q(1,8),L=2))
        obj.state=object()
        with self.assertRaises(ValueError):p.plan(0)
    def test_invalid_profile_inputs(self):
        for args in ((-1,1,2,'f0'),(0,0,2,'f0'),(0,2,1,'f0'),(0,1,2,'f1')):
            with self.assertRaises(ValueError):q112.power_recipe(*args)
        with self.assertRaises(ValueError):q112.cutoff_recipe(0,1,0)
        with self.assertRaises(TypeError):q112.KernelProfiles(None)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
