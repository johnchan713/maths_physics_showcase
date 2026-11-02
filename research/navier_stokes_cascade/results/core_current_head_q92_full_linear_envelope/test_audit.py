import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q92
class Tests(unittest.TestCase):
    def test_variable_graph_second_derivative(self):
        x,y=s.symbols('x y');v=1+x+y;f=x*y+y**3
        D=lambda z:s.diff(z,x)+v*s.diff(z,y)
        exact=D(D(f))
        frozen=s.diff(f,x,2)+2*v*s.diff(f,x,y)+v*v*s.diff(f,y,2)
        correction=(s.diff(v,x)+v*s.diff(v,y))*s.diff(f,y)
        self.assertEqual(s.expand(exact-frozen-correction),0)
        self.assertNotEqual(s.expand(correction),0)
    def test_mixed_graph_majorant(self):
        x,y=s.symbols('x y');v=1+x+y;f=(1+x+2*y)**5
        D=lambda z:s.diff(z,x)+v*s.diff(z,y)
        def norm(z,n):return sum(abs(s.diff(z,x,i,y,j-i).subs({x:0,y:0}))/s.factorial(i)/s.factorial(j-i) for j in range(n+1) for i in range(j+1))
        for j in range(4):
            N=j+2;g=1+norm(v,N)
            self.assertLessEqual(norm(D(D(f)),j),2*N*N*g*g*norm(f,N))
    def test_gain_and_ceil_frequency(self):
        for k in (Q(1,100000),Q(1,8)):
            self.assertEqual(min(q92.gains(k).values()),Q(1,2)-k)
            self.assertGreater(Q(1,2)-k,0)
        for k in (0,Q(1,4)):
            with self.assertRaises(ValueError):q92.gains(k)
        for q in range(1,15):
            eps=Q(1,q*q);self.assertLessEqual(s.ceiling(1/s.sqrt(eps)),2/s.sqrt(eps))
    def test_full_group_polynomial_and_omission_controls(self):
        data={k:Q(1) for k in q92.KEYS};out=q92.envelope(0,data)
        # d=2: E=6+2+3+15+16+10+4=56; H=4+1+4=9.
        self.assertEqual(out['F'].record(),q92.EP(56+4+9+56).record())
        self.assertEqual(set(out['groups']),set(q92.gains()))
        self.assertTrue(out['pressure'].record());self.assertTrue(out['principal_on_curl'].record())
        no_pi=q92.envelope(0,dict(data,Pi=0))
        self.assertNotEqual(out['F'].record(),no_pi['F'].record())
        self.assertEqual(out['input_order'],2);self.assertEqual(out['root_or_target_order'],3)
        self.assertEqual(out['F'].record(),q92.envelope(0,data,-1)['F'].record())
        with self.assertRaises(ValueError):q92.envelope(0,{})
        with self.assertRaises(ValueError):q92.envelope(0,data,0)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
