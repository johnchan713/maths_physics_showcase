import unittest,json
import sympy as s
import audit,q81
class Tests(unittest.TestCase):
    def test_same_column_target_cancellation(self):
        C=s.Matrix([[2,1],[1,3]]);T=s.Matrix([5,7]);h=s.Matrix([1,2,3])
        z=C.inv()*T;P=z[0]*(h*h.T)
        self.assertEqual(P,(C.inv()*T)[0]*(h*h.T))
        changed=C.inv()*(T+s.Matrix([1,0]))
        self.assertNotEqual(changed[0]*(h*h.T),P)
        self.assertEqual(changed-z,C.inv()*s.Matrix([1,0]))
    def test_moving_boundary_and_nonlocal_support(self):
        u,t=s.symbols('u t',positive=True);e=1+t
        f=u**4*(e-u)**4
        integral=s.integrate(f,(u,0,e))
        # Finite-order model only; actual support edge is smooth to all orders.
        for j in range(3):
            self.assertEqual(s.diff(f,t,j).subs(u,e),0)
        self.assertEqual(s.simplify(s.diff(integral,t)-s.integrate(s.diff(f,t),(u,0,e))),0)
        self.assertEqual(s.simplify(s.integrate(f,(u,0,e))/integral),1)
        self.assertGreater(integral.subs(t,1),0)
    def test_interval_intersection(self):
        self.assertEqual([q81.length_in_U(2,8,e) for e in (1,2,5,8,10)],[0,0,3,6,6])
        with self.assertRaises(ValueError):q81.length_in_U(8,2,5)
    def test_evidence(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
