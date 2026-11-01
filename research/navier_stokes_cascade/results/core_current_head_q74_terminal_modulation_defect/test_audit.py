import json,math,unittest
import sympy as s
import audit,q74

class Tests(unittest.TestCase):
    def test_exact_quadratic_coefficients(self):
        t,x,w,e,u,a,b=s.symbols('t x w e u a b')
        E=e*s.exp(t*a); U=u+t*b
        rows=[x*U,w*E,w*E*U,x*(U**2-E**2/2),E**2/2]
        for f,v in zip(rows,q74.quadratic(x,w,e,u,a,b)):
            self.assertEqual(s.simplify(s.diff(f,t,2).subs(t,0)/2-v),0)

    def test_normalize_before_subtracting(self):
        M,I,J,S,P,r,e,l,w=s.symbols('M I J S P r e l w',nonzero=True)
        a,b,c,d,f=q74.target((M,I,J,S,P),r,e,l,w)
        reconstructed=(-a*r*e,-c*w*e,(-a+l*b)*w*e**2,d*r*e**2,-f*e**2)
        for p,v in zip(reconstructed,(M,I,J,S,P)):self.assertEqual(s.simplify(p-v),0)

    def test_terminal_pressure_mean(self):
        # Trapezoidal diagnostic on a smooth periodic fixture, with zero
        # primitive and first slow derivative at the two integration ends.
        for N in (16,32):
            total=0.; expected=0.
            for i in range(8192):
                y=i/8192; env=math.sin(math.pi*y)**2; E=1+y/3
                A=env*math.sin(2*math.pi*N*y)
                total+=E**2*math.expm1(2*A/N)/2
                expected+=E**2*env**2/2
            self.assertLess(abs(total*N**2/expected-1),0.015)
            self.assertGreater(total,0)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
