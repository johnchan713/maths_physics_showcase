import json,unittest
from fractions import Fraction as Q
from math import factorial
import sympy as s
import audit,q72

class Tests(unittest.TestCase):
    def test_native_slope_factor(self):
        z=s.symbols('z');h=Q(1,16);f=1/(1-z/h)
        for m in range(6):
            self.assertEqual(sum(s.diff(f,z,i+1).subs(z,0)/factorial(i) for i in range(m+1)),q72.slope(m,h))

    def test_width_cost(self):
        t=Q(1,32)
        for i in range(1,6):
            raw=162*3260**i*factorial(i)**3/t**i
            self.assertEqual(q72.switch(i,t)-q72.switch(i-1,t),raw/factorial(i))

    def test_forward_densities_and_memory(self):
        x,z=s.symbols('x z',positive=True);f=1+x*z;u=x+z;E=s.sqrt(2*x)*f
        rows=[x*u,x*s.sqrt(2*x)*E,x*s.sqrt(2*x)*E*u,x*(u*u-E**2/2),E**2/2]
        major=q72.densities(x,f,u)
        for i in (0,1,2,4):self.assertEqual(s.simplify(rows[i]-major[i]),0)
        self.assertEqual(s.simplify(major[3]-rows[3]-2*x*x*f**2),0)
        # Two profiles with identical starting field and trace at x=1.
        F=f+(x-1)*z
        change=s.integrate(F**2-f**2,(x,1,x))
        self.assertEqual(change.subs(x,1),0)
        self.assertNotEqual(s.diff(change,z).subs({x:2,z:1}),0)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
