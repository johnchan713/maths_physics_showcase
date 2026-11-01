import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q86

def fixture():
    f=q86.full;r=f.refresh;p=r.maps.update.profile.pipeline;P=p.roots.P
    root=p.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
    hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
    target=r.maps.update.profile.build(root,lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3),Rmin=1,Rmax=8,
      fixed=lambda n,pole:{k:1 for k in p.fresh.FIXED if k!='rho'},compact=lambda n,pole:dict(KT=1,KA=1,I=6),slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},hull=hull,fast_hull=(1,8),max_output=14,profile=[10000]*28,profile_hull=hull)
    state=r.maps.update.SignedCovariance(target,11);current=r.maps.NewPoleMaps(state)
    return f.build(current,lambda n,p:dict(order=n,pole=p,Ez=7),0,d=1,c0=2,c1=4,shape=[1]*15,time_inverse=lambda n,p:1)
class Tests(unittest.TestCase):
    def test_actual_factory_identity_and_manufactured_leaf_trace(self):
        post=fixture();out=q86.FinalMeans(post);calls=[];EP=q86.full.EP
        # Manufactured coefficient leaves avoid evaluating the giant norm programs.
        def values(name,n,k):
            calls.append((name,n));return {v:EP(k+i) for i,v in enumerate(('x','u','z'))}
        post.current.old_full=lambda n:values('old_full',n,1)
        post.current.fast_unweighted=lambda n:values('complete_fast',n,4)
        post.current.slow=lambda n,weighted:values('fresh_slow',n,7) if weighted is False else None
        v=out.fields(0)
        self.assertEqual(calls,[('old_full',0),('complete_fast',0),('fresh_slow',0)])
        for i,k in enumerate(('x','u','z')):self.assertEqual(v[k].record(),EP(243*(12+3*i)).record())
        self.assertIs(out.state,post.current.current.state)
        with self.assertRaises(ValueError):out.fields(1)
        with self.assertRaises(TypeError):q86.FinalMeans(lambda n:1)
        post.current.fast.rows.current=object()
        with self.assertRaises(ValueError):q86.FinalMeans(post)
    def test_wrong_rate_shortcut_and_reserve(self):
        k=Q(1,100000);B=1-k;C0=2-2*k
        self.assertGreater(C0,B);self.assertGreater(C0,2-3*k)
        y=s.symbols('y',real=True);v=s.cos(2*s.pi*y)
        self.assertNotEqual(s.diff(v,y),0)
        self.assertEqual(q86.ledger(200002)['profile_covariance'],200029)
        self.assertEqual(q86.ledger(200002)['incoming_rows'],200012)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
