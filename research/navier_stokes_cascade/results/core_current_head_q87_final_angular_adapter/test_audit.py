import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q87

def fixture(M):
    f=q87.full;r=f.refresh;p=r.maps.update.profile.pipeline;P=p.roots.P
    root=p.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
    hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
    target=r.maps.update.profile.build(root,lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3),Rmin=1,Rmax=8,
      fixed=lambda n,pole:{k:1 for k in p.fresh.FIXED if k!='rho'},compact=lambda n,pole:dict(KT=1,KA=1,I=6),slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},hull=hull,fast_hull=(1,8),max_output=M+14,profile=[10000]*(M+28),profile_hull=hull)
    state=r.maps.update.SignedCovariance(target,M+11)
    return f.build(r.maps.NewPoleMaps(state),lambda n,p:dict(order=n,pole=p,Ez=7),M,d=1,c0=2,c1=4,shape=[1]*(M+15),time_inverse=lambda n,p:1)
class Tests(unittest.TestCase):
    def test_difference_matches_independent_polynomial_and_rates(self):
        old=dict(x=Q(2),u=Q(3),z=Q(5));inc=dict(x=Q(7),u=Q(11),z=Q(13))
        fixed=dict(b=Q(17),V=Q(19),G=Q(23),graph=Q(29),r=Q(31))
        d=6;rt=17*11+19*7+2*11+7*3+7*11
        at=23*11+19*13+5*11+13*3+13*11
        expected=d*11+(59*d+62)*rt+d*at+(d*d*(59**2+1)+31*59*d+31**2)*11
        self.assertEqual(q87.angular_difference(2,old,inc,fixed),expected)
        k=Q(1,100000);self.assertGreater(3-6*k,2-2*k)
    def test_chart_conversion(self):
        u,R=s.symbols('u R');f=R**3
        NR=sum(abs(s.diff(f,R,k).subs(R,4))/s.factorial(k) for k in range(4))
        NU=sum(abs(s.diff(f.subs(R,u*u),u,k).subs(u,2))/s.factorial(k) for k in range(4))
        self.assertLessEqual(NU,q87.native_factor(3,8,1)*NR)
        with self.assertRaises(ValueError):q87.native_factor(2,8,0)
    def test_exact_fast_cancellation(self):
        y=s.symbols('y',real=True);incoming=3+s.cos(y)+2*s.sin(2*y)
        fast_time=-(incoming-3);rf=7*s.cos(3*y);rs=1+s.sin(y)
        final=s.expand(incoming+fast_time+rf+rs)
        self.assertEqual(s.simplify(final-(3+rf+rs)),0)
        self.assertNotEqual(s.simplify(final-(incoming+rf+rs)),0)
    def test_typed_factory_trace_and_reserve(self):
        post=fixture(3);f=q87.means.FinalMeans(post);a=q87.FinalAngular(f);calls=[];EP=q87.full.EP
        def vector(name,n):calls.append((name,n));return {k:EP(i+1) for i,k in enumerate(('x','u','z'))}
        post.current.old_full=lambda n:vector('old',n)
        post.current.fast_unweighted=lambda n:vector('fast',n)
        post.current.slow=lambda n,w:vector('slow',n)
        post.current.initial.budget._data=lambda n,p:{k:EP(1) for k in ('b','V','G','graph','r')}
        f.state.target.source.haar.rows=lambda n:dict(Etheta=EP(1))
        # Isolate coefficient leaves too: computing a real signed norm here
        # would recurse into the deliberately replaced manufactured fixed data.
        f.state.coefficients=lambda n:{k:EP(1) for k in ('A','B','L','S')}
        f.state.covariance=lambda n:dict(complete=EP(7),pole=20)
        out=a.coefficient(0)
        self.assertEqual(calls,[('old',3),('fast',3),('slow',3)])
        self.assertEqual(out['Etheta'].record(),(243*(out['Haar']+out['fast']+out['slow'])).record())
        self.assertEqual(out['covariance_request'],1)
        with self.assertRaises(ValueError):a.coefficient(1)
        with self.assertRaises(TypeError):q87.FinalAngular(lambda n:1)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
