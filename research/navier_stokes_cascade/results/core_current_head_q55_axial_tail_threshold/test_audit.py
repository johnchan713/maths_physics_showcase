import json
import unittest
from fractions import Fraction as F
import sympy as s
import audit
import threshold as t


class Tests(unittest.TestCase):
    def test_geometric_sums_between_dyadic_shells(self):
        for m in (4,5,199989,199990):
            for r in (F(1),F(17,3),F(64),F(127,2)):
                hi,lo=t.split(m)
                self.assertGreater(2*hi-m,0)
                self.assertLess(2*lo-m,0)
                self.assertLess(sum((2**h/r)**(2*(hi if 2**h<=r else lo)-m)
                                    for h in range(25)),4)

    def test_pell_divisor_and_shift_resonance(self):
        U,a,y,z,M=s.symbols('U a y z M')
        for j in (2,5,12,29):
            A,B=t.pell(j);x=A-B; w=x+(1-s.sqrt(2))*B
            self.assertEqual(s.expand(w*(A+s.sqrt(2)*B)),(-1)**j)
            phase=lambda u,y,z:x*y+B*z-M*w*u
            self.assertEqual(s.expand(phase(U+a,y+M*a,z+(1-s.sqrt(2))*M*a)-phase(U,y,z)),0)
            self.assertNotEqual(s.diff(phase(U+a,y,z),a),0)

    def test_integer_threshold(self):
        k=F(1,100000)
        for m in (199989,199990):
            rate=2+m*k/2; before=rate-k/2
            if m==199989:
                self.assertLessEqual(before,3-6*k);self.assertGreater(rate,3-6*k)
            else:
                self.assertLess(before,3-5*k);self.assertEqual(rate,3-5*k)

    def test_replay(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
