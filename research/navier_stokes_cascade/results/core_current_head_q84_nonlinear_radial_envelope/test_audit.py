import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q84
class Tests(unittest.TestCase):
    def test_exact_projected_remainder(self):
        y=s.symbols('y',nonzero=True);ep=s.symbols('ep');b,V,x,u,z,a,a1,ts,dz,vis=s.symbols('b V x u z a a1 ts dz vis')
        co=(y+y**-1)/2;si=(y-y**-1)/(2*s.I)
        beta=x*si;v=u*co;gamma=z*si
        qc=lambda f:s.expand(f).expand().collect(y)-s.expand(f).coeff(y,0)
        def inv(f):
            f=s.expand(qc(f));return s.expand(sum(f.coeff(y,k)*y**k/(s.I*k) for k in range(-6,7) if k))
        G,rr,tt,rt,rz,tz,prr,ptt,prt,E=s.symbols('G rr tt rt rz tz prr ptt prt E')
        Wrr=rr*co;Wtt=tt*co;Wrt=rt*si;Wrz=rz*co;Wtz=tz*si
        Prr=prr*co;Ptt=ptt*co;Prt=prt*si
        J=b*v+V*beta+beta*v+Wrt
        F=-E*co-v+dz*(G*v+V*gamma+gamma*v+Wtz)-vis*v
        H=b*gamma+G*beta+beta*gamma+Wrz
        rest=-ts*beta-(2*b*beta+beta**2+Wrr)-dz*H+Wtt+ep*dz**2*beta
        original=qc(-a1*inv(J)-a*inv(F)+v**2+rest)
        Drr=Wrr-ep*Prr/2;Dtt=Wtt-ep*Ptt/2;Drt=Wrt-ep*Prt/2
        rem=-a1*inv(b*v+V*beta+beta*v+Drt)-a*inv(F)+v**2-ts*beta-(2*b*beta+beta**2)-dz*H+Dtt-Drr+ep*dz**2*beta
        principal=ep*(qc(Ptt-Prr)-a1*inv(Prt))/2
        self.assertEqual(s.expand(original-principal-qc(rem)),0)
        self.assertNotEqual(s.expand(qc(Dtt-Drr)),0)
    def test_rates_and_request_ceiling(self):
        k=Q(1,100000);r=q84.rates(k)
        self.assertEqual(r['F'],2-3*k);self.assertGreaterEqual(r['nonF'],2-2*k)
        data={k:Q(1) for k in ('b0','V','G','x','u','z','W','Z','e','r','g','c','T')}
        out=q84.envelope(199996,data)
        self.assertEqual(out['field_order'],200002);self.assertEqual(out['angular_row_order'],200000)
        recipe=q84.angular.temporal_recipe(199996)
        self.assertEqual(recipe['raw_source_order'],200000);self.assertFalse(recipe['materialized'])
        more=dict(data,Z=2);self.assertGreater(q84.envelope(0,more)['projected_remainder'],q84.envelope(0,data)['projected_remainder'])
        with self.assertRaises(ValueError):q84.envelope(0,dict(data,Z=-1))
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
