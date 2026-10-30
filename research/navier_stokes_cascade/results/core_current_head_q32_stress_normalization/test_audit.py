import unittest,json
from fractions import Fraction as Q
from math import factorial
import sympy as s
import normalization as n
import audit

class NormalizationTests(unittest.TestCase):
    def test_physical_covariance_exponents(self):
        h=s.symbols('h');A=s.Rational(1,2)+h;alpha=1+h
        self.assertEqual(s.expand(2*A-alpha),h)
        self.assertEqual(s.expand(-2*A+h+alpha),0)
        self.assertNotEqual(s.expand(-A+h+alpha),0)
        self.assertEqual(s.expand(-2*A+3*h),h-1)

    def test_missing_prefactor_has_nonzero_derivatives(self):
        t=s.symbols('t',positive=True);h=s.Rational(1,1000)
        for j in range(1,5):
            self.assertNotEqual(s.diff(t**(-1-h),t,j).subs(t,1),0)

    def test_independent_slow_derivative_fixture(self):
        R,Z,T=s.symbols('R Z T');scale=T+Z**2
        # h=0 coordinate limiting fixture only; no assertion h=0 is admissible.
        f=(1+R**2/(2*scale)+(Z/s.sqrt(scale))**2)/scale
        for m in range(4):
            for i in range(m+1):
                for j in range(m+1-i):
                    k=m-i-j
                    value=abs(s.diff(f,R,i,Z,j,T,k).subs({R:2,Z:0,T:1}))
                    self.assertLessEqual(value,8*n.physical_cost(m,1,2))

    def test_radial_majorant_against_differentiation(self):
        U=s.symbols('U',positive=True)
        for d in (Q(1,2),Q(1),Q(2),Q(3)):
            f=(U**(1/s.Rational(d)))**3
            for m in range(5):
                v=abs(s.diff(f,U,m).subs(U,1))
                self.assertLessEqual(v,48*n.radial_cost(m,d,1,2))
        self.assertEqual(n.radial_cost(0,2,1,2),1)
        with self.assertRaises(ValueError):n.radial_cost(2,2,0,2)

    def test_factorial_count(self):
        for m in range(6):
            total=Q(0)
            for a in range(m+1):
                for b in range(m+1-a):
                    for c in range(m+1-a-b):
                        for d in range(m+1-a-b-c):
                            e=m-a-b-c-d
                            total+=Q(1,factorial(a)*factorial(b)*factorial(c)*factorial(d)*factorial(e))
            self.assertEqual(total,Q(5**m,factorial(m)))

    def test_adapter_into_original_root_recurrence(self):
        P=n.roots.P
        target=lambda j:n.transfer(j,n.stress.fixture(j),1,2)['root_target']
        r=n.roots.RootBudget(H=lambda j:P(1),target=target,wave=lambda j:P(1),
                            curl=lambda j:P(1),loghalf=lambda j:P(1),
                            K0=1,beta=1,R=1,mass=2,width=Q(1,4))
        self.assertEqual(r.B(0).evaluate(1),target(0).evaluate(1))
        self.assertEqual(r.B(1).evaluate(1),(target(1)+target(0)).evaluate(1))
        self.assertEqual(r.data('target',2).evaluate(1),target(2).evaluate(1))

    def test_higher_stress_order_and_retained_pole(self):
        b=n.stress.fixture(2);v=n.transfer(2,b,1,2)
        self.assertEqual(v['higher_stress_raw'],4*n.physical_cost(2,1,2)*b['correction'])
        self.assertEqual((v['higher_epsilon_order'],v['full_epsilon_order'],v['pole']),(3,1,9))
        with self.assertRaises(ValueError):n.transfer(1,b,1,2)

    def test_overlap_derivatives_cancel_only_after_sum(self):
        x=s.symbols('x',real=True);f=s.exp(x)
        terms=[s.sin(x)**2*f,s.cos(x)**2*f]
        self.assertNotEqual(s.diff(terms[0],x),s.sin(x)**2*s.diff(f,x))
        for j in range(4):self.assertEqual(s.simplify(sum(s.diff(t,x,j) for t in terms)-s.diff(f,x,j)),0)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['all_source_premises_closed'])

if __name__=='__main__':unittest.main()
