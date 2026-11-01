import unittest,json
import sympy as s
import audit,q77

class Tests(unittest.TestCase):
    def test_affine_end_energy_root(self):
        x=s.symbols('x');p=x+2*x*x;A=s.Integer(1)
        f=lambda a:a-a*a/8-(a+p)**2/8-(s.Rational(3,4)+x+x*x)
        for _ in range(3):
            A=(A-f(A)/(1-A/4-(A+p)/4)).series(x,0,5).removeO().expand()
        self.assertEqual(q77.energy_root(4,1,2,s.Rational(1,8),[0,1,1,0,0],[(s.Rational(1,8),1,[0,1,2])]),[A.coeff(x,j) for j in range(5)])

    def test_heat_energy_units(self):
        l,Xp,Xe,Ep,Ee,delta=s.symbols('l Xp Xe Ep Ee delta',positive=True)
        # Incoming S changes by delta, future integral E² by 2 delta.
        ratio=Xe*Ee**2/(Xp*Ep**2)
        self.assertEqual(s.simplify(l*delta/(Xp*Ep**2)-l*ratio*(2*delta/(Xe*Ee**2))/2),0)

    def test_relative_angular_shape(self):
        z=s.symbols('z',real=True);theta=s.Rational(2,3); f=1/(1+z*z)
        normalized=(f**theta/f)**2
        self.assertEqual(s.simplify(normalized-(1+z*z)**(2*(1-theta))),0)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
