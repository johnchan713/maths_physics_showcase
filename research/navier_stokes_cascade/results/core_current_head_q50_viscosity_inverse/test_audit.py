import json
import unittest
from fractions import Fraction as Q
import sympy as s
import inverse as v
import audit


class Tests(unittest.TestCase):
    def test_symbol_and_pell_witness(self):
        b=s.sqrt(2)-1
        for n,(x,y) in enumerate(v.pell_modes(12),1):
            l=x-b*y;t=b*x+y
            self.assertEqual(s.expand(l*t-(-1)**n*b),0)
            self.assertTrue(bool(t>=(abs(x)+abs(y))/s.Integer(3)))
            self.assertEqual(s.simplify((l**3/t)*t**4-(-1)**n*b**3),0)

    def test_inverse_on_real_modes(self):
        import numpy as np
        # Exact zero at the origin, two real cosines plus a constant.
        y,z=np.meshgrid(np.arange(32)/32,np.arange(32)/32,indexing='ij')
        b=2**0.5-1
        f=2-np.cos(2*np.pi*(y-2*z))-np.cos(2*np.pi*(3*y+z))
        self.assertEqual(f[0,0],0)
        h=np.fft.fftn(f); m=np.fft.fftfreq(32)*32
        x,w=np.meshgrid(m,m,indexing='ij')
        l=x-b*w;t=b*x+w; active=(abs(h)>1e-8)&((x!=0)|(w!=0))
        g=np.zeros_like(h);g[active]=-(2*np.pi)**2*l[active]**3/t[active]*h[active]
        recovered=np.zeros_like(h);recovered[active]=-g[active]*t[active]/((2*np.pi)**2*l[active]**3)
        self.assertAlmostEqual(-np.fft.ifftn(recovered)[0,0].real,f.mean())
        self.assertLess(np.max(abs(np.fft.ifftn(recovered).real-(f-f.mean()))),1e-12)

    def test_physical_pullback(self):
        k=Q(1,100000); eta=Q(1,200000)
        self.assertEqual(v.target_order(3-6*k+eta),2-3*k+eta)
        for rho in (Q(1,6),Q(1,3),Q(2,3),Q(5,6)):
            c=Q(1,64); j=Q(3,7)
            sharp=v.INVERSE_CONSTANT*Q(2)**int(6*(3*rho-1))/j**3
            self.assertLessEqual(sharp,v.necessary_target(1,j,1/c))
        with self.assertRaises(ValueError):v.necessary_target(1,0,4)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))


if __name__=='__main__':unittest.main()
