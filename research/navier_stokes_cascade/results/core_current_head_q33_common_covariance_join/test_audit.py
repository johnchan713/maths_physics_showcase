import unittest,json
from fractions import Fraction as Q
import sympy as s
import covariance as c
import audit

class CovarianceTests(unittest.TestCase):
    def test_profile_coordinates_and_weight_invariant(self):
        r,R,Z,scale=s.symbols('r R Z scale',positive=True);D=s.Rational(499,1000)
        self.assertEqual(s.simplify((s.sqrt(r)*R)**2/(2*r*scale)-R**2/(2*scale)),0)
        self.assertEqual(s.simplify(r**D*Z*(r*scale)**(-D)-Z*scale**(-D)),0)

    def test_full_real_covariance_not_just_selected_fluxes(self):
        t=s.symbols('t',real=True);A=s.Matrix([1,2,3]);B=s.Matrix([4,-2,1]);k=3
        v=A*s.cos(t)-B*s.sin(t)/k
        W=(A*A.T+B*B.T/k**2)/2
        for i in range(3):
            for j in range(i,3):
                avg=s.integrate(v[i]*v[j],(t,0,2*s.pi))/(2*s.pi)
                self.assertEqual(avg,W[i,j])
                self.assertNotEqual(W[i,j],A[i]*A[j]/2)

    def test_norm_adapter_preserves_exact_polynomials_and_adds_stress_once(self):
        p=c.normal.roots.fixture();g=dict(Rmin=1,Rmax=2)
        for n in range(3):
            v=c.common_inputs(n,p,c.normal.stress.fixture,g)
            old=c.haar.weighted_covariance(n,c.PrimaryProxy(p))
            self.assertEqual(v['W'].c,old['W'].c)
            self.assertEqual(v['D'].c,(old['Dcurl']+v['higher']).c)
            H=c.normal.transfer(n,c.normal.stress.fixture(n),**g)['higher_stress_raw']
            self.assertEqual(v['higher'].evaluate(1),H)
            self.assertEqual(v['pole'],12*n+18)

    def test_exact_unchanged_row_interface(self):
        v=c.common_inputs(1,c.normal.roots.fixture(),c.normal.stress.fixture,dict(Rmin=1,Rmax=2))
        data=c.haar.fixture();data.update({k:v[k] for k in ('W','D')})
        rows=c.haar.row_budget(0,data)
        self.assertEqual(rows['field_order'],3)
        self.assertGreater(rows['Ez'].evaluate(1),0)

    def test_curl_request_and_scale_rescaling(self):
        trace=set();p=c.normal.roots.fixture(trace)
        c.common_inputs(2,p,c.normal.stress.fixture,dict(Rmin=1,Rmax=2))
        self.assertIn(('target',3),trace)
        self.assertNotIn(('target',4),trace)
        proxy=c.PrimaryProxy(p)
        self.assertEqual(proxy.amplitude(1).c,p.amplitude(1).c)

    def test_orders_and_poles(self):
        for n in range(8):self.assertLessEqual(3*n+3,12*n+18)
        for k in (Q(1,100000),Q(1,8)):
            for e in (Q(1),Q(1,2),Q(1,100)):
                # Exact comparison after raising to a positive common denominator.
                exponent=2-2*k
                self.assertLessEqual(e**(3*exponent.denominator),e**exponent.numerator)
        with self.assertRaises(ValueError):c.common_inputs(0,None,None,{},kappa=0)

    def test_original_canonical_factory_preserves_value_data(self):
        normal=c.columns.entrance.reference.normal
        p=normal.frame.inputs(bmin=Q(1,20),bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),cmin=1,cmax=1,Dbox=6,Mz=1)
        sp=normal.selector.fixture();L=normal.normal_length_certificate(p)['length_rule_floor']
        sp['length']=Q(-(-L.numerator//L.denominator))
        ell=normal.selector.select_tail(sp)['ell_tail']
        P=c.columns.P
        kw=dict(F=lambda j:P(j+2),G=lambda j:P(j+3),btilde=lambda j:P(j+1),
                path=lambda j:P({1:1}) if j==0 else P(int(j==1)),tilt=Q(1,100))
        _,column,_,_,_=c.columns.from_original(p,sp,ell,**kw)
        template=c.columns.entrance.inverse.roots.fixture()
        joined=c.canonical_root(template,column,lambda j:c.columns.P(j+1),
                                c.normal.stress.fixture,dict(Rmin=1,Rmax=2))
        self.assertEqual((joined.K0,joined.beta,joined.width),(template.K0,template.beta,template.width))
        self.assertIs(joined.curl,template.curl)
        expected=c.normal.transfer(0,c.normal.stress.fixture(0),1,2)['root_target']
        self.assertEqual(joined.target(0).c,expected.c)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['full_source_or_residual_or_blowup'])

if __name__=='__main__':unittest.main()
