from fractions import Fraction as Q
from itertools import product
import unittest,json
import cone as c
import audit

class ConeTests(unittest.TestCase):
    def test_cramer_identity_and_bounds(self):
        p=(Q(2),Q(0));m=(Q(0),Q(2));t=(Q(3,5),Q(4,5))
        r=c.enclosure(p,m,t,(Q(1,100),)*3);b=c.root_constants([r],2)
        for signs in product((-1,1),repeat=6):
            v=[tuple(x[j]+Q(signs[2*i+j],200) for j in range(2)) for i,x in enumerate((p,m,t))]
            P,M,T=v;D=c.det(P,M);a=c.det(T,M);z=c.det(P,T)
            roots=(2*a/D,2*z/D)
            self.assertEqual(tuple(P[j]*roots[0]+M[j]*roots[1] for j in range(2)),tuple(2*x for x in T))
            self.assertGreaterEqual(min(roots),b['beta'])
            inv=max(abs(M[1])+abs(M[0]),abs(P[1])+abs(P[0]))/abs(D)
            self.assertLessEqual(inv,b['K0'])
            for actual,ref,E in zip((D,a,z),(c.det(p,m),c.det(t,m),c.det(p,t)),r['error_bounds']):
                self.assertLessEqual(abs(actual-ref),E)

    def test_orientation_must_be_common(self):
        r=c.enclosure((0,2),(2,0),(Q(3,5),Q(4,5)),(0,0,0),-1)
        self.assertGreater(c.root_constants([r],1)['beta'],0)
        rbad=c.enclosure((0,2),(2,0),(Q(3,5),Q(4,5)),(0,0,0),1)
        with self.assertRaises(ValueError):c.root_constants([rbad],1)

    def test_invertibility_does_not_supply_positive_root(self):
        for t in ((Q(3,5),Q(-4,5)),(1,0)):
            r=c.enclosure((1,0),(0,1),t,(0,0,0))
            self.assertEqual(r['margins'][0],1)
            with self.assertRaises(ValueError):c.root_constants([r],1)

    def test_positive_pulse_mass_does_not_supply_column_independence(self):
        pulse=(1,1,0);self.assertGreater(sum(x*x for x in pulse),0)
        col=(pulse[0]*pulse[1],pulse[0]*pulse[2])
        self.assertEqual(c.det(col,col),0)
        r=c.enclosure(col,col,(1,0),(0,0,0))
        with self.assertRaises(ValueError):c.root_constants([r],1)

    def test_quadratic_errors_are_required(self):
        r=c.enclosure((0,0),(0,0),(0,0),(2,3,5))
        self.assertEqual(r['error_bounds'],(6,15,10))
        self.assertEqual(c.det((2,0),(0,3)),6)

    def test_cover_aggregation(self):
        r1=c.enclosure((2,0),(0,2),(Q(3,5),Q(4,5)),(0,0,0))
        r2=c.enclosure((3,0),(0,3),(Q(4,5),Q(3,5)),(0,0,0))
        b=c.root_constants([r1,r2],Q(1,2))
        self.assertEqual(b['d'],4);self.assertEqual(b['a'],Q(6,5));self.assertEqual(b['M'],3)
        self.assertEqual(b['beta'],Q(1,15));self.assertEqual(b['K0'],Q(3,2))

    def test_unbounded_rho_preserves_lower_bound(self):
        r=c.enclosure((1,0),(0,1),(Q(3,5),Q(4,5)),(0,0,0));b=c.root_constants([r],1)
        for n in (1,2,10000):self.assertGreaterEqual(n**3*Q(3,5),b['beta'])

    def test_guards(self):
        with self.assertRaises(TypeError):c.enclosure((1.,0),(0,1),(1,1),(0,0,0))
        with self.assertRaises(ValueError):c.enclosure((1,0),(0,1),(1,1),(-1,0,0))
        with self.assertRaises(ValueError):c.root_constants([],1)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['actual_head_margins_proved'])

if __name__=='__main__':unittest.main()
