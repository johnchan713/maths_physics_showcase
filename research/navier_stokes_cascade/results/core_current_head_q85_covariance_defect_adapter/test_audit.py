import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q85

def fixture(M):
    r=q85.full.refresh;p=r.maps.update.profile.pipeline;P=p.roots.P
    root=p.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
    hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
    t=r.maps.update.profile.build(root,lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3),Rmin=1,Rmax=8,
      fixed=lambda n,pole:{k:1 for k in p.fresh.FIXED if k!='rho'},compact=lambda n,pole:dict(KT=1,KA=1,I=6),
      slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},hull=hull,fast_hull=(1,8),max_output=M+3,profile=[10000]*(M+17),profile_hull=hull)
    return r.maps.update.SignedCovariance(t,M)
def coeff(x):return {j:dict(v.terms) for j,v in x.c.items()}
class Tests(unittest.TestCase):
    def test_exact_real_tensor_expansion(self):
        ep,q=s.symbols('epsilon q');A=s.Matrix([1,2,3]);B=s.Matrix([2,-1,4]);L=s.Matrix([-1,3,2]);S=s.Matrix([3,1,-2])
        sym=lambda x,y:x*y.T+y*x.T
        W=ep*(A+ep*L)*(A+ep*L).T/2+ep**2*q**2*(B+ep*S)*(B+ep*S).T/2
        defect=ep**2*sym(A,L)/2+ep**3*L*L.T/2+ep**2*q**2*B*B.T/2+ep**3*q**2*sym(B,S)/2+ep**4*q**2*S*S.T/2
        self.assertEqual((W-ep*A*A.T/2-defect).applyfunc(s.expand),s.zeros(3))
        self.assertNotEqual((defect-ep**3*L*L.T/2).applyfunc(s.expand),defect.applyfunc(s.expand))
        self.assertNotEqual(A*L.T,L*A.T)
    def test_original_state_adapter_and_rate_guard(self):
        state=fixture(0);root=state.target.root;out=q85.defect(state,0,target_rate=Q(2))
        self.assertIs(state.primary.primary.root,root)
        self.assertEqual(coeff(out['Z']),coeff(out['oldcurl']+243*state.covariance(0)['complete']))
        self.assertEqual(out['pole'],state.covariance(0)['pole']);self.assertEqual(out['target_order_required'],3)
        self.assertEqual(out['raw_to_factorial'],243)
        self.assertTrue(coeff(out['oldcurl']));self.assertTrue(coeff(out['signed']))
        with self.assertRaises(ValueError):q85.defect(state,0,target_rate=2-Q(3,100000))
        with self.assertRaises(ValueError):q85.defect(state,1,target_rate=Q(2))
        with self.assertRaises(TypeError):q85.defect(object(),0,target_rate=Q(2))
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
