from fractions import Fraction as Q
import unittest,json
import audit
import lift as t
class LiftTests(unittest.TestCase):
    def test_path_time_and_transverse(self):
        for b in (Q(1,10),Q(2,5),Q(9,10)):
            H=(Q(3),Q(-2));c=Q(1,256);scale=Q(17);offset=Q(7);radius=Q(1,19)
            for w in (0,1,2*radius/c):
                p=t.path(H,w,b,c,scale,offset,radius)
                self.assertEqual(p[0]-b*p[1],H[0]-b*H[1])
                self.assertEqual((scale*(b*p[0]+p[1])/(1+b*b)-offset+radius)/c,w)
    def test_path_composition(self):
        H=(2,3);b=Q(2,5);c=Q(1,300)
        self.assertEqual(t.path(t.path(H,7,b,c),11,b,c),t.path(H,11,b,c))
    def test_exact_affine_derivatives(self):
        b=Q(2,5);H=(Q(2),Q(3));c=Q(1,400);scale=3;h=Q(1,7)
        base=t.path(H,5,b,c,scale);J=t.jacobian(b)
        for j in (0,1):
            X=list(H);X[j]+=h;p=t.path(X,5,b,c,scale)
            self.assertEqual(tuple((p[i]-base[i])/h for i in (0,1)),tuple(J[i][j] for i in (0,1)))
        p=t.path(H,5+h,b,c,scale)
        self.assertEqual(tuple((p[i]-base[i])/h for i in (0,1)),(c*b/scale,c/scale))
    def test_projection_and_column_bounds(self):
        for b in (Q(1,100),Q(2,5),Q(99,100)):
            J=t.jacobian(b)
            for i in (0,1):
                self.assertEqual(J[i][0]*b+J[i][1],0)
                for j in (0,1):self.assertEqual(sum(J[i][k]*J[k][j] for k in (0,1)),J[i][j])
            for j in (0,1):self.assertLessEqual(sum(abs(J[i][j]) for i in (0,1)),2)
    def test_shifted_lift_not_averaged(self):
        a=t.path((0,0),1,Q(2,5),Q(1,256),offset=0)
        b=t.path((0,0),1,Q(2,5),Q(1,256),offset=1)
        self.assertNotEqual(a,b)
    def test_jet_factor(self):
        self.assertEqual(t.jet_factor(3,2,256),Q(32,256**2))
        self.assertEqual(t.jet_factor(7,0,10**30),128)
    def test_guards(self):
        with self.assertRaises(ValueError):t.jet_factor(-1,0,1)
        with self.assertRaises(ValueError):t.jet_factor(0,0,0)
        with self.assertRaises(ValueError):t.path((1,2),0,0,1)
    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['actual_extended_source_bounds_proved'])
if __name__=='__main__':unittest.main()
