from fractions import Fraction as Q
from functools import lru_cache
import json
import unittest
import sympy as s
import waves as w
import audit

class Phase:
    def normal(self,j):return w.P(2 if j==0 else 1 if j==1 else 0)
    def matrix(self,j):return w.P()
    def moving_normal(self,j):return w.P(int(j==0))

def fixture(trace=None,xi_trace=None,cut_slope=1):
    root=w.roots.RootBudget(H=lambda j:w.P(int(j==0)),target=lambda j:w.P(int(j==0)),
        wave=lambda j:w.P(int(j==0)),curl=lambda j:w.P(int(j==0)),
        loghalf=lambda j:w.P(),K0=1,beta=1,R=1,mass=1,width=Q(1,8),trace=trace)
    def xi(j):
        if xi_trace is not None:xi_trace.add(j)
        return w.P(int(j==0))
    return w.SignedWaveInputs(root,cut_wave=lambda j:w.P(1 if j==0 else cut_slope if j==1 else 0),
        xi=xi,xi_pole=3,phase=Phase(),inverse_normal_square=1)

class WaveTests(unittest.TestCase):
    def test_same_root_all_value_and_nonwave_slots(self):
        a=fixture();r,c=a.regular,a.cut
        for name in ('H','target','curl','loghalf','K0','beta','R','mass','width'):
            self.assertEqual(getattr(r,name),getattr(c,name))
        self.assertIsNot(r.wave,c.wave)
        self.assertEqual(r.amplitude(0).c,c.amplitude(0).c)

    def test_cutoff_derivative_is_in_actual_U_only(self):
        a,b=fixture(cut_slope=0).through(0),fixture(cut_slope=3).through(0)
        for name in ('T','R','Pi'):self.assertEqual(a[name].c,b[name].c)
        self.assertGreater(b['U'].evaluate(1),a['U'].evaluate(1))

    def test_independent_cutoff_curl_identity(self):
        x,y,z=s.symbols('x y z');A=s.Matrix([0,0,x*x*y]);psi=1+y
        def curl(v):return s.Matrix([s.diff(v[2],y)-s.diff(v[1],z),s.diff(v[0],z)-s.diff(v[2],x),s.diff(v[1],x)-s.diff(v[0],y)])
        comm=s.simplify(curl(psi*A)-psi*curl(A))
        self.assertEqual(comm,s.Matrix([x*x*y,0,0]))
        self.assertNotEqual(comm,s.zeros(3,1))

    def test_moving_normal_pressure_and_reciprocal_jets(self):
        a=fixture();x=s.symbols('x');row=-1/(1+x*x)
        # n=(x,1,0), K=0, nprime=(1,0,0), |n|^2>=1 on [0,1].
        for j in range(7):
            exact=abs(s.diff(row,x,j).subs(x,0))
            self.assertLessEqual(exact,a.pressure_multiplier(j).evaluate(1))
        self.assertGreater(a.through(0)['Pi'].evaluate(1),0)
        class Dropped(Phase):
            def moving_normal(self,j):return w.P()
        a.phase=Dropped()
        self.assertEqual(a.pressure_multiplier(0).evaluate(1),0)

    def test_pressure_product_derivative_cannot_be_frozen(self):
        x=s.symbols('x');t=2+x;row=1/(1+x*x)
        self.assertNotEqual(s.diff(row*t,x,2).subs(x,0), (row*s.diff(t,x,2)).subs(x,0))
        self.assertEqual(s.diff(row*t,x,2).subs(x,0),-4)

    def test_signed_square_survives(self):
        theta=s.symbols('theta',real=True);old=2*s.cos(theta);new=3*s.cos(theta)
        change=s.integrate((old+new)**2-old**2,(theta,0,2*s.pi))/(2*s.pi)
        cross=s.integrate(2*old*new,(theta,0,2*s.pi))/(2*s.pi)
        self.assertEqual(change,Q(21,2));self.assertEqual(change-cross,Q(9,2))

    def test_exact_values_and_two_sign_mass(self):
        a=fixture(cut_slope=0)
        # At S=1: primary=1, signed=1/2, curl=3/2 for each sign.
        out=a.through(0)
        self.assertEqual({k:v.evaluate(1) for k,v in out.items()},{'T':3,'R':3,'Pi':3,'U':6})
        for k,v in a.source_wave(0).items():self.assertEqual(v.c,out[k].c)

    def test_exponents_and_poles(self):
        for k in (Q(1,100000),Q(1,10),Q(249,1000)):
            alpha=2-3*k;mu=1-k
            self.assertGreaterEqual(alpha-Q(1,2),Q(1,2))
            self.assertGreaterEqual(alpha-k,mu);self.assertGreaterEqual(mu,Q(1,2))
        a=fixture()
        for n in range(7):
            for j in range(n+1):
                for p in (6*j+3,6*j+9,6*j+a.xi_pole,6*j+a.xi_pole+6):
                    self.assertGreaterEqual(a.pole(n),p)

    def test_unmodified_prefix_factory_and_full_request_depth(self):
        trace=set();xi=set();a=fixture(trace,xi)
        base=w.prefix.fixture();wave_requests=set()
        def wave(n):wave_requests.add(n);return a.source_wave(n)
        base.wave=wave
        source=base.raw_source(0)
        self.assertGreater(source.evaluate(1),0)
        self.assertEqual(max(wave_requests),15)
        self.assertEqual(max(xi),16)
        self.assertEqual(max(j for name,j in trace if name=='wave'),16)
        self.assertGreater(base.source(0)['complete_wave_wave'].evaluate(1),0)
        self.assertGreater(base.source(0)['complete_mean_wave'].evaluate(1),0)
        self.assertGreater(base.rows(1,True)['pm'].evaluate(1),0)

    def test_guards(self):
        a=fixture()
        with self.assertRaises(ValueError):a.through(-1)
        with self.assertRaises(ValueError):a.pole(True)
        with self.assertRaises(ValueError):w.SignedWaveInputs(a.regular,cut_wave=a.cut.wave,xi=a.xi,xi_pole=-1,phase=Phase(),inverse_normal_square=1)
        with self.assertRaises(ValueError):w.SignedWaveInputs(a.regular,cut_wave=a.cut.wave,xi=a.xi,xi_pole=1,phase=Phase(),inverse_normal_square=0)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['complete_first_source_or_residual_proved'])
if __name__=='__main__':unittest.main()
