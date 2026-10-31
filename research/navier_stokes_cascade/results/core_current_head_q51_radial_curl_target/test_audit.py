from fractions import Fraction as Q
import json, unittest
import target_rate as q
import audit

class Tests(unittest.TestCase):
    def test_cylindrical_split(self):
        import sympy as s
        R,Z,y,M,eps=s.symbols('R Z y M eps')
        h=s.Matrix([R,Z+y,R*y]);A=s.Matrix([Z*y,R+Z,y**2])
        C=h.cross(A);er=s.Matrix([1,0,0]);ez=s.Matrix([0,0,1])
        B=er.cross(C.diff(R)+M*R*C.diff(y))+eps*ez.cross(C.diff(Z))+ez*C[1]/R
        self.assertEqual(s.expand(B[0]+eps*(h.diff(Z).cross(A)+h.cross(A.diff(Z)))[1]),0)
        self.assertNotEqual(s.expand(B[0]),0)
        self.assertNotEqual(s.expand(B[1].diff(M)),0)
        self.assertNotEqual(s.expand(B[2].diff(M)),0)
    def test_same_polynomial_and_guard(self):
        class Lift:
            def lift(self,j):return q.curl.P({0:j+1,1:1})
        op=q.curl.CurlBudget(Lift(),1,8)
        p=q.profile;P=p.pipeline.roots.P
        for j in range(6):
            upper=op.normalized(j,'U').c
            self.assertTrue(all(upper.get(i,0)>=v for i,v in q.radial_operator(op,j).c.items()))
        def target(operator):
            one=lambda j:P(int(j==0))
            root=p.pipeline.roots.RootBudget(H=one,target=one,wave=one,curl=operator,
                loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
            hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
            return p.build(root,lambda n:dict(target=1,correction=1,pole=3*n+3),
                Rmin=1,Rmax=8,fixed=lambda n,v:{k:1 for k in p.pipeline.fresh.FIXED if k!='rho'},
                compact=lambda n,v:dict(KT=1,KA=1,I=6),slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},
                hull=hull,fast_hull=(1,8),max_output=0,profile=[10000]*14,profile_hull=hull)
        # Reserve the checked radial operator only through order 3; later fixture
        # jets are manufactured lower norm inputs, not physical curl bounds.
        t=target(lambda j:P(op.normalized(j,'U').c) if j<=3 else P(1))
        for e in (1,2):
            v=q.target_bound(t,op,0,e);self.assertIs(v['normalized_target'],t.native(0,e))
            self.assertTrue(v['normalized_target'].c)
            self.assertEqual((v['epsilon_order'],v['pole']),(2,t.pole))
        with self.assertRaises(ValueError):q.target_bound(target(lambda j:P(1)),op,0,1)
    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
