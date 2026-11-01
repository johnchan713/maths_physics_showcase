import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q79

class Tests(unittest.TestCase):
    def test_exact_ambient_lift_and_normal_memory(self):
        v=s.symbols('v',real=True);n=s.Matrix([1+v,2,3-v]);m=n.diff(v)
        K=s.Matrix([[0,-2,0],[5,0,0],[7,0,0]])
        N=(n.T*n)[0];P=n*n.T/N;I=s.eye(3);h=n/N
        A0=-(I-P)*K-n*m.T/N
        p,d=s.symbols('p d',real=True);A=A0-d*(I-P)+p*P
        zero=lambda M:self.assertEqual(M.applyfunc(s.simplify),s.zeros(*M.shape))
        zero(m.T+n.T*A-p*n.T)
        zero(h.diff(v)-(I-2*P)*m/N)
        self.assertEqual(s.simplify((h.diff(v).T*h.diff(v))[0]-(m.T*m)[0]/N**2),0)
        zero(A0*h-h.diff(v)+(I-P)*(K*n+m)/N)
        # Without moving-normal pressure the constraint fails identically.
        self.assertNotEqual(s.simplify((m.T+n.T*(A+n*m.T/N)-p*n.T)[0]),0)

    def test_energy_bound_keeps_damping_signed(self):
        # Independent exact PSD check, for a moving normal and nonsymmetric K.
        n=s.Matrix([2,1,2]);m=s.Matrix([1,-2,1]);N=(n.T*n)[0]
        P=n*n.T/N;I=s.eye(3);K=s.Matrix([[0,-2,0],[5,0,0],[7,0,0]])
        # Frobenius norm K < 9, |m| < 3, |n|=3; Dref=2.
        bound=s.Rational(12)
        for d in (0,1,10000):
            lam=s.Rational(1,3);dref=2;p=lam-dref
            A=-(I-P)*K-n*m.T/N-d*(I-P)+p*P
            B=bound*I-((A+A.T)/2-p*I)
            for idx in ((0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)):
                self.assertGreaterEqual(B.extract(idx,idx).det(),0)
        head=q79.head_budget(9,1,2,64,3,2,1,1,2,1)
        self.assertEqual(head['log_ambient'],64*(21+18*18+Q(13,2)))

    def test_finite_head_rounding_and_fixed_low_requests(self):
        # Exact scale fixture: h=1/200, ell0=800, epsilon^-1/2<4.
        kh=1+2**2
        for k in (1,2,3,4):
            for j in (-3,-1,1,3):
                self.assertGreaterEqual(abs(Q(j,k))/2,Q(1,kh*2))
        trace=set();bg=q79.background.fixture(trace);low=q79.low_inputs(bg,Q(1,10))
        self.assertEqual(trace,{('leading',1),('leading',2),('cutoff',2),('profile',1,1)})
        self.assertEqual(low['Db'],3*Q(1,10)*low['leading_two']+low['error_one'])
        self.assertGreater(bg.bound(1)['head']['correction'],bg.bound(1)['head']['tail'])
        with self.assertRaises(ValueError):q79.head_budget(0,1,2,1,1,1,1,1,1,1)
        with self.assertRaises(TypeError):q79.low_inputs(bg,0.1)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
