from fractions import Fraction as Q
import unittest,json
import directed as d
import audit

class Tests(unittest.TestCase):
    def setUp(self):self.s=d.Shell(1,8,2,4,Q(1,2))
    def test_collar(self):
        for p in (0,1,30,100003):
            r=self.s.collar(p)
            self.assertLessEqual(p*r*r,self.s.a)
            self.assertLessEqual(r,self.s.gap)
            self.assertLessEqual(4*r,self.s.lo)
    def test_weight_and_directed_ratios(self):
        import mpmath as m
        f=lambda q:m.mpf(q.numerator)/q.denominator
        with m.workdps(50):
            L=m.log(8);a=f(self.s.a);b=4
            for p in (0,3,12):
                w=lambda x:m.exp(-a/x**2-b/(L-x)**2)*min(1,x,L-x)**(-p) if 0<x<L else 0
                r=f(self.s.collar(p))
                for x in (r/2,r,L/2,L-r,L-r/2):
                    self.assertLessEqual(w(x),m.exp(p*p/(4*a)))
                    if x<=r:self.assertLessEqual(m.quad(w,[0,x])/w(x),r*m.exp(b/(f(self.s.lo)-r)**2))
                    if x>=L-r:self.assertLessEqual(m.quad(w,[x,L])/w(x),r*m.exp(a/(f(self.s.lo)-r)**2))
    def test_adapter(self):
        def jets(n,p):return dict(order=n,pole=p,cutoff_order=n+1,bounds=dict(chi=2,one_minus_chi=2,chi_U=3))
        def bad(n,p):return dict(jets(n,p),bounds=dict(kernel=1))
        with self.assertRaises(ValueError):self.s.remaining(bad)(0,2)
        v=d.radial.build(2,1,8,self.s.remaining(jets),lambda n,p:1,2,3)(1,3)
        self.assertEqual((v['order'],v['pole']),(1,3))
        self.assertGreater(len(v['KT'].c[0].terms),1)
        self.assertGreater(max(v['KA'].c[0].terms),0)
    def test_evidence(self):self.assertEqual(audit.generate(),json.loads((d.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
