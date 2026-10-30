import unittest,json
from fractions import Fraction as Q
import sympy as s
import compact as c
import audit

def data():return dict(f=2,g=3,k=5,chi=1,one_minus_chi=1,chi_U=7,jacobian=11,length=13,edge=17,kernel=19,interior=23)

class CompactTests(unittest.TestCase):
    def test_constants_independent_substitution(self):
        d=data();r=c.constants(2,**d)
        self.assertEqual(r['KT'],64*3*2*2*11*19)
        self.assertEqual(r['KA'],64*5*7*2*11*13*17*23)

    def test_logarithmic_jacobian(self):
        x,z=s.symbols('x z',real=True);U=8*z*z*s.exp(3*x/2)
        self.assertEqual(s.diff(U,x),s.Rational(3,2)*U)
        # A constant source integrates against J dx, not dx alone.
        self.assertEqual(s.simplify(s.integrate(s.diff(U,x),(x,0,2))-8*z*z*(s.exp(3)-1)),0)

    def test_shift_fixed_before_differentiation(self):
        U,H,M,a,Z=s.symbols('U H M a Z')
        shifted=(U+a)**3*(H+M*a)+Z*(U+a)
        self.assertEqual(s.diff(shifted,U),3*(U+a)**2*(H+M*a)+Z)
        self.assertEqual(s.diff(shifted,Z),U+a)
        self.assertEqual(s.diff(shifted,H), (U+a)**3)
        self.assertNotEqual(s.diff(shifted,a),s.diff(shifted,U))

    def test_directed_primitive_requires_other_edge_exclusion(self):
        x=s.symbols('x');f=x*x*(1-x)**2
        left=s.integrate(f,(x,0,1));self.assertEqual(left,s.Rational(1,30))
        # Even a vanishing source leaves nonzero left primitive at the right edge.
        self.assertNotEqual(left,0);self.assertEqual(f.subs(x,1),0)

    def test_floor_on_scaled_support_not_union(self):
        # R/sqrt(s) in [2,3] stays inside shell (1,4); union R is [1,6].
        for root_s in (Q(1,2),1,2):
            for profile_r in (2,3):
                R=root_s*profile_r
                self.assertTrue(1<R/root_s<4)
        self.assertGreater(Q(6)/Q(1,2),4)

    def test_collar_and_archive(self):
        self.assertEqual(c.collar(2,Q(1,10)),Q(1,10))
        rows,total=c.archive([6,12,6,16,32],lambda p:{'KT':Q(100,p),'KA':Q(p)})
        self.assertEqual(list(rows),[6,12,16,32])
        self.assertGreater(total['KT'],rows[32]['KT'])
        self.assertEqual(total['KA'],66)

    def test_length_is_not_double_counted(self):
        d=data();a=c.constants(0,**d);d['length']*=2;b=c.constants(0,**d)
        self.assertEqual(a['KT'],b['KT']);self.assertEqual(2*a['KA'],b['KA'])

    def test_guards(self):
        with self.assertRaises(ValueError):c.collar(1,1)
        with self.assertRaises(ValueError):c.constants(-1,**data())
        with self.assertRaises(ValueError):c.archive([],lambda p:None)
        d=data();d['jacobian']=1.0
        with self.assertRaises(TypeError):c.constants(0,**d)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['slow_response_and_full_source_closed'])
if __name__=='__main__':unittest.main()
