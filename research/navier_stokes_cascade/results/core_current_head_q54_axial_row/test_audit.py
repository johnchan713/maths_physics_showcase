import json
import unittest
from fractions import Fraction as F
import sympy as s
import audit
import axial


class Axial(unittest.TestCase):
    def test_two_viscosities_and_complement(self):
        # A nonzero Fourier mode: graph and time have different symbols.
        # Keep a free tail jet: its value and its viscous image both survive.
        R,M,e,t=s.symbols('R M e t',positive=True)
        j=R+1; l=2*s.I; I=1/(3*s.I); W=R**2+1
        graph=lambda f:s.diff(f,R)+M*j*l*f
        visc=lambda f:-e*(graph(graph(f))+graph(f)/R)
        f0=graph(W)+W/R
        f1=visc(-I*f0)
        tail=t*(R**2+R)
        z1=-I*(f1-tail)
        after=f1+z1/I+visc(z1)
        poly=s.Poly(after.expand(),M)
        self.assertEqual(s.simplify(poly.LC()-e**2*j**5*l**5*I**2*W),0)
        self.assertEqual(poly.degree(),5)
        self.assertEqual(s.simplify(s.diff(after,t)*t-tail-visc(I*tail)),0)
        self.assertNotEqual(s.expand(visc(I*tail)),0)

    def test_rate_separation_and_least_count(self):
        k=F(1,100000);d=axial.rates(k)
        self.assertEqual(d['leading']-d['required'],k)
        self.assertEqual(d['remainder']-d['leading'],k)
        self.assertLess(d['direct_tail'],d['required'])
        self.assertLess(2+(d['p']-1)*k,d['leading'])
        self.assertEqual(2+d['p']*k,d['leading'])
        self.assertEqual(d['phi_extra'],199992)

    def test_replay(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
