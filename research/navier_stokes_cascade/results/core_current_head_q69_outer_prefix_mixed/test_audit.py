import json,unittest
import sympy as s
from math import factorial
import audit,q69

class Tests(unittest.TestCase):
    def test_local_log_taper_operator(self):
        t=s.symbols('t', positive=True); u=s.symbols('u');f=u**5-2*u**3
        for n in range(1,6):
            arg=s.log(1+t)/64
            recovered=sum((-1)**(n-j)*v*s.diff(f,u,j).subs(u,arg)/64**j for j,v in enumerate(q69.stirling(n)))/(1+t)**n
            self.assertEqual(s.simplify(s.diff(f.subs(u,arg),t,n)-recovered),0)

    def test_exponential_recurrence(self):
        t=s.symbols('t')
        H=[0,1,3,0,2,1,0]
        actual=s.exp(sum(h*t**j for j,h in enumerate(H))).series(t,0,7).removeO()
        self.assertEqual(q69.exponential(H),[actual.coeff(t,j) for j in range(7)])

    def test_physical_change_and_taper(self):
        X,y=s.symbols('X y',positive=True)
        for i in range(1,7):
            f=sum((-1)**j*y**j for j in range(i+1))
            actual=s.diff(f.subs(y,s.log(X/2)),X,i).subs(X,2)/factorial(i)
            # Alternating derivatives attain the unsigned Stirling sum.
            self.assertEqual(abs(actual),q69.physical(i,2,[factorial(j) for j in range(i+1)]))
        self.assertEqual(q69.taper(1),s.Rational(q69.step(1),16))
        self.assertEqual(q69.taper(2),4*(s.Rational(q69.step(1),64)+s.Rational(q69.step(2),64**2)))

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
