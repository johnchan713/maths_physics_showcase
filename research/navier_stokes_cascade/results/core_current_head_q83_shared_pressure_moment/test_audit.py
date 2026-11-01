import unittest,json
import sympy as s
import audit,q83
class Tests(unittest.TestCase):
    def test_independent_radial_adjoint_and_axial_weights(self):
        u,z=s.symbols('u z',positive=True);R=u**3;D=R.diff(u)/R
        F=u*z+z*z;G=u*u*z;V=R*F;C=1+z*z
        rr,tt,zz,rt,rz=s.symbols('rr tt zz rt rz')
        a1=6*V/R**2-2*V.diff(u)/(R.diff(u)*R)
        J=D*(tt-rr)-R.diff(u)*a1*rt
        K=D*R**2*(zz-(rr+tt)/2)+2*R*G.diff(u)*rz-R**2*F.diff(u)*rt
        A=D*R**2*(zz-rr)+D*(C-R**2/2)*(tt-rr)+2*R*G.diff(u)*rz+((2*C-R**2)*F.diff(u)-4*C*D*F)*rt
        self.assertEqual(s.expand(A-K-C*J),0)
        self.assertEqual(s.expand(A.diff(z)-K.diff(z)-C*J.diff(z)),s.expand(C.diff(z)*J))
        self.assertNotEqual(s.expand(C.diff(z)*J),0)
    def test_request_from_both_shell_branches(self):
        for m,z,p in [(3,1,99995),(4,0,199995)]:
            for r in (0,1):
                rows=[(m+z+k,2*p-r-k) for k in (p+1,p-1)]
                got=q83.request(m,z,p,r)
                self.assertEqual(got['total'],max(a+b for a,b in rows))
                self.assertEqual(got['profile'],1+max(a for a,b in rows))
        self.assertEqual(q83.record()['radial_extra_output']['eta'],2)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
