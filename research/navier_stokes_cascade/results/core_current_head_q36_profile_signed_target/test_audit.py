import unittest,json
from fractions import Fraction as Q
from math import factorial
import profile_target as p
import audit

def fixture(trace=None):
    P=p.pipeline.roots.P
    root=p.pipeline.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),
        wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),
        K0=1,beta=1,R=1,mass=2,width=Q(1,4))
    stress=lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3)
    hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
    kw=dict(Rmin=1,Rmax=8,fixed=lambda N,pole:{k:1 for k in p.pipeline.fresh.FIXED if k!='rho'},
        compact=lambda n,pole:dict(KT=1,KA=1,I=6),
        slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},
        hull=hull,fast_hull=(1,8),max_output=1,trace=trace,
        profile=[10000]*15,profile_hull=hull)
    return root,stress,kw

class ProfileTargetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.trace=set();root,stress,kw=fixture(cls.trace)
        cls.v=p.build(root,stress,**kw)
        cls.source0=cls.v.source.source(0,1)
    def test_projector_mass_and_moving_pressure(self):
        import sympy as s
        x=s.symbols('x');f=30*(x-2)**2*(3-x)**2
        self.assertEqual(s.integrate(f,(x,2,3)),1)
        a=s.integrate(f*f,(x,2,3));moment=s.integrate(x*x*f,(x,2,3))/2
        self.assertGreater(a,0)
        for scale in (Q(1,4),Q(1),Q(4)):
            # R=sqrt(s)x makes R^e sigma_e dR exactly f^2/A2 dx.
            for e in (1,2):
                R=s.sqrt(scale)*x
                sigma=R**(-e)/s.sqrt(scale)*f*f/a
                self.assertEqual(s.integrate(R**e*sigma*s.sqrt(scale),(x,2,3)),1)
            self.assertEqual(s.integrate((s.sqrt(scale)*x)**2*f/2,(x,2,3)),scale*moment)
        self.assertNotEqual(moment,0) # c_rho,Z P cannot be omitted away from Z=0.
    def test_independent_profile_and_projector_derivatives(self):
        import sympy as s
        R,T=s.symbols('R T',positive=True);x=R/s.sqrt(T)
        f=30*(x-2)**2*(3-x)**2;rho=f/s.sqrt(T)
        z=s.symbols('z');A=s.integrate((30*(z-2)**2*(3-z)**2)**2,(z,2,3))
        c=self.v.source.profile
        for n in range(4):
            def norm(expr):
                return sum(abs(s.diff(expr,R,i,T,j).subs({R:s.Rational(5,2),T:1}))/(factorial(i)*factorial(j))
                           for i in range(n+1) for j in range(n+1-i))
            self.assertLessEqual(norm(rho),c.rho(n))
            for e in (1,2):self.assertLessEqual(norm(R**(-e)*s.sqrt(T)*rho**2/A),c.sigma(n,e))
        self.assertGreater(c.dilation(2),0)
        self.assertGreater(c.c_rho(1),c.c_rho(0))
    def test_exact_symbolic_ring_transfer(self):
        f=p.pipeline.fresh
        v=f.EP({2:f.ExpSum([(Q(2,3),Q(1000000)),(7,Q(1,5))]),4:f.ExpSum([(9,0)])})
        w=p.rebind(v,p.signed.projector.haar.sharp.fresh)
        self.assertEqual(v.record(),w.record())
        self.assertEqual(p.rebind(w,f).record(),v.record())
    def test_one_profile_density_pressure_projector(self):
        v=self.v;c=v.source.profile;h=v.source.haar
        self.assertIs(h.c_rho.__self__,c)
        self.assertEqual(h.initial.budget.fixed(3,h.pole)['rho'].record(),c.weighted_rho(3).record())
        n=0;t=h.initial.output(1)
        expected=t['Jz']+c.c_rho(1)*t['P']
        self.assertEqual(self.source0['moment'].record(),expected.record())
        angular=v.source.source(n,2)
        self.assertEqual(angular['moment'].record(),t['Jtheta'].record())
        self.assertEqual(v.source.bump(0,1).record(),(c.floor*v.source.pull.pullback(0,1)*c.sigma(0,1)).record())
    def test_original_root_target_and_full_curl(self):
        v=self.v
        self.assertIs(v.root,v.source.haar.initial.budget.primary.primary.root)
        self.assertEqual(v.ratio,8)
        self.assertEqual(v.pole,2*(12*14+18))
        self.assertTrue(v.signed_amplitude(0).record())
        expected=v.signed_amplitude(0)+v.signed_amplitude(1)
        self.assertEqual(v.signed_curl(0).record(),expected.record())
        self.assertEqual(max(n for k,n in self.trace if k=='covariance'),14)
        with self.assertRaises(ValueError):v.signed_curl(1)
        self.assertFalse(v.ledger()['actual_signed_wave_full_covariance_or_residual_certified'])
    def test_input_boundaries(self):
        root,stress,kw=fixture()
        for change in ({'profile':[1]*14},{'c_rho':lambda n:1},
                       {'profile_hull':dict(kw['profile_hull'],t1=Q(1,2))}):
            with self.assertRaises(ValueError):p.build(root,stress,**dict(kw,**change))
        bad=p.build(root,stress,**dict(kw,fixed=lambda n,pole:{k:1 for k in p.pipeline.fresh.FIXED}))
        with self.assertRaises(ValueError):bad.source.source(0,1)
        with self.assertRaises(ValueError):self.v.source.profile.sigma(0,3)
    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
