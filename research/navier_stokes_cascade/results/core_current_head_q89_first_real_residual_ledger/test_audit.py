import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q89
class Tests(unittest.TestCase):
    def test_projection_orthogonality_and_mean_pressure(self):
        th,R,z=s.symbols('theta R z',real=True);avg=lambda f:s.integrate(f,(th,0,2*s.pi))/(2*s.pi)
        v=s.Matrix([2+s.cos(th),3+s.sin(2*th),R+z+s.cos(3*th)])
        zero=v.applyfunc(avg);non=v-zero
        self.assertEqual(s.simplify(avg(v.dot(v))-zero.dot(zero)-avg(non.dot(non))),0)
        p=R**2*z+z**3;grad=s.Matrix([s.diff(p,R),0,s.diff(p,z)])
        self.assertEqual(grad-grad.applyfunc(avg),s.zeros(3,1))
        self.assertGreater(avg(non.dot(non)),0)
    def test_spatial_pressure_cutoff_is_not_temporal_tail(self):
        R,th,t=s.symbols('R theta t');p=s.cos(th)
        self.assertEqual(s.diff(R*p,R)-R*s.diff(p,R),p)
        self.assertEqual(s.diff(t*p,R)-t*s.diff(p,R),0)
    def test_stage_and_rate_guards(self):
        with self.assertRaises(ValueError):q89.check_ledger('first_sourced',dict(q89.SOURCES),[])
        wrong=dict(q89.SOURCES,axial_compact='final_radial_first_real')
        with self.assertRaises(ValueError):q89.check_ledger('first_real',wrong,[])
        with self.assertRaises(ValueError):q89.check_ledger('first_real',dict(q89.SOURCES),['spatial'])
        for h in (Q(1,200),Q(1,1000)):
            r=q89.rates(h,m=2);k=Q(1,100000)
            self.assertEqual(r['missing_source_gain'],2-5*k)
            self.assertEqual(r['physical_total'],-Q(5,2)-h*(1+k)-4)
            self.assertLess(r['normalized_total'],r['target'])
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
