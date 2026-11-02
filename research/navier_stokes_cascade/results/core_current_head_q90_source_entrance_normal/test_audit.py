import unittest,json
from math import comb
import sympy as s
import audit,q90
class Tests(unittest.TestCase):
    def test_moving_normal_complete_pressure(self):
        v=s.symbols('v',real=True);n=s.Matrix([1,v,0]);t=s.Matrix([-v*(1+v),1+v,v**2])
        K=s.Matrix([[0,2,1],[3,0,v],[1,2,0]]);pi=1+v**3;omega=2;d=1+v
        f=-t.diff(v)-K*t-d*t-s.I*omega*n*pi
        recovered=s.I*((K.T*n-n.diff(v)).dot(t)+n.dot(f))/(omega*n.dot(n))
        self.assertEqual(s.simplify(n.dot(t)),0)
        self.assertEqual(s.simplify(recovered-pi),0)
        omitted=s.I*((K.T*n).dot(t)+n.dot(f))/(omega*n.dot(n))
        self.assertNotEqual(s.simplify(omitted-pi),0)
    def test_exact_entrance_recurrence_and_one_order_shift(self):
        v=s.symbols('v');H=1+v+v**2
        for r in range(5):
            t=v**(r+1)*(1+v);g=H*t-s.diff(t,v)
            J=lambda f:[s.diff(f,v,j).subs(v,0) for j in range(6)]
            out=q90.entrance_jets(J(H),J(g))
            self.assertEqual(out[:7],[s.diff(t,v,j).subs(v,0) for j in range(7)])
            self.assertEqual(out[:r+1],[0]*(r+1))
            self.assertNotEqual(out[r+1],0)
    def test_normal_only_forcing_cannot_be_dropped(self):
        n=s.Matrix([1,0,0]);f=n*3;P=s.eye(3)-n*n.T
        self.assertEqual(P*f,s.zeros(3,1))
        self.assertEqual(s.I*n.dot(f)/2,3*s.I/2)
        self.assertGreater(q90.pressure([0],[1],[0],[3],2)[0],0)
    def test_flat_exit_source_does_not_give_zero_exit_velocity(self):
        v=s.symbols('v');f=v**2*(1-v)**2
        for x in (0,1):
            self.assertEqual(f.subs(v,x),0);self.assertEqual(s.diff(f,v).subs(v,x),0)
        self.assertEqual(-s.integrate(f,(v,0,1)),-s.Rational(1,30))
    def test_mixed_pressure_envelope(self):
        v,x=s.symbols('v x');a=s.Matrix([1+v+x,(1+v+x)**2,0]);b=s.Matrix([1,v,x]);u=s.Matrix([1+x,v**2,1]);g=s.Matrix([1+v,x**2,v+x])
        def bound(z,j):
            return max(sum(abs(s.diff(y,v,i,x,j-i).subs({v:0,x:0})) for y in z) for i in range(j+1))
        arrays=[[bound(z,j) for j in range(5)] for z in (a,b,u,g)]
        out=q90.pressure(*arrays,2);p=(a.dot(u)+b.dot(g))/2
        for j in range(5):
            for i in range(j+1):self.assertLessEqual(abs(s.diff(p,v,i,x,j-i).subs({v:0,x:0})),out[j])
        self.assertEqual(q90.pressure([0],[1],[0],[3],-2),[s.Rational(3,2)])
        with self.assertRaises(ValueError):q90.pressure([0],[1],[0],[3],0)
        with self.assertRaises(ValueError):q90.product([1],[1,2])
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
