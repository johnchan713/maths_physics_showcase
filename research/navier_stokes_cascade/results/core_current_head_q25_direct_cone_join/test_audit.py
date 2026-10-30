from fractions import Fraction as Q
from itertools import product
import unittest,json
import direct as d
import audit

def det(p,m):return p[0]*m[1]-p[1]*m[0]
class DirectTests(unittest.TestCase):
    def test_all_error_corners_and_cone_endpoints(self):
        for g in (Q(1,100),Q(1,2),Q(1)):
            e=g/4
            for a,b,c,f in product((-e,e),repeat=4):
                p=(1+a,1+b);m=(1+c,-1+f);D=-det(p,m)
                self.assertGreaterEqual(D,2-4*e-2*e*e)
                self.assertLessEqual(D,2+4*e+2*e*e)
                for r in (-(1-g),1-g):
                    self.assertGreaterEqual(-det((1,r),m),g-e*(2-g))
                    self.assertGreaterEqual(-det(p,(1,r)),g-e*(2-g))

    def test_mass_cancellation_and_original_root_floor(self):
        for jp,jm in ((Q(1,10),Q(10)),(Q(10),Q(1,10))):
            g=Q(1,2);b=d.bounds(g,g/4,2,3,Q(1,10),10,2)
            for signs in product((-1,1),repeat=4):
                a,z,c,f=(g*x/4 for x in signs)
                p=(jp*2*(1+a),jp*3*(1+z));m=(jm*2*(1+c),jm*3*(-1+f))
                t=(Q(4),Q(3)) # Vhat=2,r=1/2
                D=det(p,m);roots=(det(t,m)/D,det(p,t)/D)
                self.assertGreaterEqual(min(roots),b['beta'])
            self.assertGreaterEqual(b['beta'],b['original_beta'])

    def test_original_inverse_and_beta_constants(self):
        for g in (Q(1,10000),Q(1,2),Q(1)):
            b=d.bounds(g,g/4,2,3,Q(1,10),10,2)
            self.assertLess(b['inverse_bound'],b['original_K0'])
            self.assertGreaterEqual(b['beta'],b['original_beta'])

    def test_orthogonal_basis_cost_is_not_one(self):
        O=((Q(3,5),Q(-4,5)),(Q(4,5),Q(3,5)))
        for i in (0,1):
            for j in (0,1):self.assertEqual(sum(O[k][i]*O[k][j] for k in (0,1)),int(i==j))
        norm=max(sum(abs(x) for x in r) for r in O)
        self.assertGreater(norm,1);self.assertLess(norm,2)

    def test_universal_quarter_bounds(self):
        e=Q(1,4)
        self.assertEqual(2-4*e-2*e*e,Q(7,8));self.assertEqual(2+4*e+2*e*e,Q(25,8))
        self.assertLess(Q(40,7),6);self.assertGreater(Q(4,25),Q(1,8))

    def test_selector_implies_direction_gate(self):
        for C,amin,u,g,r0 in ((Q(7),Q(1,3),Q(2),Q(1,4),Q(1,6057072)),):
            length=(4*C*max(1/amin,1/u)/g)**2
            # Integer upper bound avoids materializing any physical small q.
            n=int(length/(2*r0))+1
            self.assertGreaterEqual(2*r0*n*n,length)
            self.assertLessEqual(C*C*max(1/amin,1/u)**2/length,g*g/16)

    def test_guards(self):
        with self.assertRaises(ValueError):d.bounds(Q(1,2),Q(1,4),1,1,1,1,1)
        with self.assertRaises(TypeError):d.bounds(.5,0,1,1,1,1,1)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['historical_minimum_band_membership_proved'])
if __name__=='__main__':unittest.main()
