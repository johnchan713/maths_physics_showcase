"""Exact model controls and guarded same-chart constant substitution."""
from fractions import Fraction as Q
import unittest,json
import sympy as s
import uniform as u
import audit

class UniformTests(unittest.TestCase):
    def test_homogeneous_and_forced_Gaussian_equations(self):
        v,L=s.symbols('v L',positive=True)
        P=s.exp(-(v-L/2)**2/(2*L));d=s.Rational(3,2)+v/L
        generator=s.diag(3-d,-3-d)
        for y,force in ((s.exp(v),0),(s.exp(v)-1,1)):
            z=s.Matrix([P*y,0]);rhs=generator*z+s.Matrix([P*force,0])
            self.assertEqual((s.diff(z,v)-rhs).applyfunc(s.simplify),s.zeros(2,1))
        self.assertEqual((P*(s.exp(v)-1)).subs(v,0),0)
        self.assertEqual(P.subs(v,L/2),1)
        self.assertEqual(s.simplify(s.diff(P,v)/P-(2-d)),0)

    def test_arbitrary_reference_envelope_identity(self):
        v=s.symbols('v');P=s.Function('P')(v);lam=s.Function('lambda')(v)
        dref=lam-s.diff(P,v)/P;z=P*s.exp(v)
        self.assertEqual(s.simplify(s.diff(z,v)-(lam+1-dref)*z),0)

    def test_propagator_on_all_subintervals(self):
        v,w,t=s.symbols('v w t');V=lambda a,b:s.diag(s.exp(a-b),s.exp(-5*(a-b)))
        self.assertEqual((V(v,t)*V(t,w)-V(v,w)).applyfunc(s.simplify),s.zeros(2))
        self.assertEqual(s.diff(V(v,w),v),s.diag(1,-5)*V(v,w))
        self.assertEqual(V(w,w),s.eye(2))

    def test_ordinary_coefficients_are_uniform_but_scaled_error_is_not(self):
        v,L=s.symbols('v L',positive=True);d=s.Rational(3,2)+v/L
        self.assertEqual(s.diff(d,v),1/L);self.assertEqual(s.diff(d,v,2),0)
        for n in (1,4,100):
            self.assertEqual(d.subs({v:0,L:n}),s.Rational(3,2))
            self.assertEqual(d.subs({v:n,L:n}),s.Rational(5,2))
        self.assertEqual(s.integrate(1,(v,0,L)),L)

    def test_fixed_reserve_cannot_absorb_model_loss(self):
        for h,delta,M in ((Q(1,200),3,7),(Q(1,10**100),10**20,9),(Q(1,3),0,0)):
            c=h*delta+2*M;ell=max(Q(1),2*c+1)
            self.assertGreaterEqual(u.reserve_log_lower(ell,h,delta,M),ell**2/2)
            self.assertGreater(u.reserve_log_lower(2*ell,h,delta,M),u.reserve_log_lower(ell,h,delta,M))

    def test_three_coordinate_mean_value_factor_is_needed(self):
        R,Z,T=s.symbols('R Z T');F=R+Z+T+R*(Z+T)
        # F_R varies in two slow directions; F itself varies in all three.
        linear=R+Z+T
        for a in (Q(1,100),Q(1,7)):
            self.assertEqual(linear.subs({R:a,Z:a,T:a}),3*a)
        self.assertEqual(s.diff(F,R),1+Z+T)

    def test_background_requests_and_exact_Db_formula(self):
        trace=set();b=u.head.fixture(trace);r=u.required(b,Q(1,8))
        self.assertEqual(r['Db'],Q(3,8)*b.leading_bound(2)+b.bound(1)['error_FG'])
        self.assertEqual(r['Mz'],b.bound(1)['full_FG'])
        self.assertIn(('profile',1,1),trace);self.assertIn(('cutoff',2),trace)
        self.assertIn(('leading',2),trace)
        self.assertFalse(any(x[0]=='profile' and x[1]>1 for x in trace))

    def test_old_selector_not_mutated_and_dominance_failure_rejected(self):
        b=u.head.fixture();p=u.frame.inputs(Rmax=2);saved=dict(p)
        r=u.join(b,p);self.assertEqual(p,saved);self.assertEqual(r['original'],saved)
        self.assertFalse(all(r['dominance'].values()))
        with self.assertRaises(ValueError):u.require_reuse(r)

    def test_guarded_recipe_uses_original_constants(self):
        b=u.head.fixture();need=u.required(b,Q(1,100));p=u.frame.inputs(Rmax=2,Dbox=Q(1,100),Db=2*need['Db'],Mz=3*need['Mz'])
        r=u.join(b,p);c=u.require_reuse(r);old=u.frame.budgets(p)
        self.assertTrue(all(r['dominance'].values()))
        self.assertEqual(c['H'],old['CE']+old['Cd'])
        self.assertGreaterEqual(c['power_of_3'],c['H'])
        self.assertLess(c['power_of_3'],c['H']+1)
        self.assertEqual(c['scale_degree'],0)
        self.assertTrue(c['original_value_gates_still_required'])
        self.assertEqual(r['original'],p)

    def test_each_dominance_slot_is_required(self):
        b=u.head.fixture();need=u.required(b,1)
        for key in ('Db','Mz'):
            inputs=dict(need);inputs[key]/=2
            r=u.join(b,u.frame.inputs(Rmax=2,Dbox=1,**inputs))
            self.assertFalse(r['dominance'][key])
            with self.assertRaises(ValueError):u.require_reuse(r)

    def test_guards(self):
        with self.assertRaises(ValueError):u.join(u.head.fixture(),u.frame.inputs(Rmax=3))
        with self.assertRaises(ValueError):u.required(u.head.fixture(),-1)
        with self.assertRaises(ValueError):u.reserve_log_lower(0,1,1,1)
        with self.assertRaises(TypeError):u.reserve_log_lower(2,0.1,1,1)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['actual_original_constant_dominance_proved'])
        self.assertFalse(r['symbolic']['actual_extended_source_bounds_proved'])

if __name__=='__main__':unittest.main()
