import importlib.util
import json
from math import comb, factorial
from fractions import Fraction as Q
import unittest
import sympy as s
import stress_inputs as f
import audit


class StressTests(unittest.TestCase):
    def test_stirling_operator_on_independent_monomials(self):
        x = s.symbols('x', positive=True)
        for j in range(1, 8):
            row = f.stirling_abs(j)
            for power in (s.Rational(3, 7), 2, 11):
                rhs = sum((-1)**(j-k)*c*power**k for k, c in enumerate(row))
                self.assertEqual(s.simplify(s.diff(x**power,x,j)*x**(j-power)),rhs)

    def test_bell_recurrence_by_independent_series(self):
        log, edge = f.weight_jets(7,Q(1,16),Q(1,2))
        z=s.symbols('z')
        p=sum(s.Rational(v)*z**j/factorial(j) for j,v in enumerate(log))
        # A finite exponential series suffices at every requested degree.
        out=s.Poly(sum(p**k/factorial(k) for k in range(8)).expand(),z)
        for j in range(8):
            self.assertEqual(out.nth(j)*factorial(j),edge[j])

    def test_actual_two_edge_weight_and_reciprocal(self):
        x=s.symbols('x', positive=True);a=s.Rational(1,16)
        log=-a/(1+s.log(x))**2-4/(1-s.log(x))**2
        _,E=f.weight_jets(5,Q(1,16),Q(1,2))
        for sign in (-1,1):
            p=s.Integer(1)
            for j in range(6):
                self.assertLessEqual(abs(p.subs(x,1)),E[j])
                p=s.diff(p,x)+sign*s.diff(log,x)*p
        # Outer constant four must contribute even if the inner width is tiny.
        self.assertEqual(f.weight_jets(1,Q(1,10**6),1)[0][1],Q(2,10**6)+8)

    def test_geometric_tail_identity_and_endpoint(self):
        t=s.symbols('t',positive=True)
        for N in range(2,8):
            exact=s.summation((t/2)**s.Symbol('n',integer=True),
                              (s.Symbol('n',integer=True),N,s.oo))
            # Independent finite sum with its exact remaining tail.
            partial=sum(Q(1,2**n)*t**(n-2) for n in range(N,N+20))
            remain=Q(1,2**(N+20))*t**(N+18)/(1-t/2)
            formula=Q(1,2**N)*t**(N-2)/(1-t/2)
            self.assertEqual(s.cancel(partial+remain-formula),0)
            self.assertEqual(formula.subs(t,1),Q(1,2**(N-1)))
            for tv in (Q(1,1000),Q(1,2),Q(99,100)):
                self.assertLess(formula.subs(t,tv),Q(1,2**(N-1)))

    def test_exceptional_first_quotient_is_not_compact(self):
        y=s.symbols('y',positive=True)
        self.assertEqual(s.limit(y**-3*s.exp(s.Rational(1,16)/(2-y)**2),y,0),s.oo)
        with self.assertRaises(ValueError):f.compact_tail(3,{1:7,2:8})
        self.assertEqual(f.compact_tail(0,{}),Q(1,2))
        self.assertEqual(f.compact_tail(4,{2:7,3:8}),Q(121,8))

    def test_mixed_cutoff_weighted_stress_fixture(self):
        x,q=s.symbols('x q',positive=True);h=Q(1,1000)
        log=-s.Rational(1,16)/(1+s.log(x))**2-4/(1-s.log(x))**2
        poly=1+x+x*x;chi=(1-q)**3
        _,edge=f.weight_jets(3,Q(1,16),Q(1,2))
        pnorm=[Q(19,4),4,2,0] # raw polynomial bounds on [1/2,3/2]
        B=[max(sum(Q(comb(j,i))*edge[i]*pnorm[j-i] for i in range(j+1))
               for j in range(m+1)) for m in range(4)]
        K=[];ce=chi
        for j in range(4):
            K.append(max([Q(1)]+[Q(sum(abs(c) for c in s.Poly(ce,q).all_coeffs()))]))
            ce=q*s.diff(ce,q)
        K=[max(K[:j+1]) for j in range(4)]
        b=f.stress_budget(3,h,Q(1,16),Q(1,2),lambda j:B[j],lambda j:B[j],
                          lambda j:10,lambda j:K[j],lambda n,j:1)
        # Differentiate after cancelling only the common value q^(2h)*zeta.
        for k in range(4):
            radial=poly
            for _ in range(k):radial=s.diff(radial,x)+s.diff(log,x)*radial
            for j in range(4-k):
                temporal=chi
                for _ in range(j):temporal=q*s.diff(temporal,q)+2*s.Rational(h)*temporal
                value=abs((radial*temporal).subs({x:1,q:s.Rational(1,2)}))
                self.assertLessEqual(value,b['correction'])

    def test_scalar_chart_has_no_velocity_prefactor(self):
        path=audit.HERE.parent/'core_background_finite_head_budget/budgets.py'
        spec=importlib.util.spec_from_file_location('q31_background',path)
        bg=importlib.util.module_from_spec(spec);spec.loader.exec_module(bg)
        self.assertEqual(bg.chart_constants(0,1,2)['composition'],1)
        R,Z,T=s.symbols('R Z T');ss=T+Z**2
        # h=0 is a coordinate-only limiting fixture, not the construction h.
        expr=(R**2/(2*ss))**2+(Z/s.sqrt(ss))**2+s.log(ss)**2
        # Raw (ell,X,eta) jets through order two at this base are <=8.
        for i in range(3):
            for j in range(3-i):
                for k in range(3-i-j):
                    value=abs(s.diff(expr,R,i,Z,j,T,k).subs({R:2,Z:0,T:1}))
                    bound=8*factorial(2)*bg.chart_constants(2,1,2)['composition']
                    self.assertLessEqual(value,bound)

    def test_positive_h_and_derivative_guards(self):
        args=dict(a=Q(1,16),xlo=1,leading=lambda j:1,first_weighted=lambda j:1,
                  first_relative=lambda j:1,cutoff=lambda j:1,compact_constants=lambda n,j:1)
        for h in (0,Q(1,100),-1):
            with self.assertRaises(ValueError):f.stress_budget(1,h=h,**args)
        with self.assertRaises(ValueError):f.stress_budget(1,h=Q(1,1000),derivative_basis='frozen',**args)
        self.assertEqual(f.stress_budget(1,h=Q(1,10**100),**args)['pole'],6)

    def test_evidence(self):
        r=audit.generate()
        self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['full_source_or_residual_proved'])


if __name__=='__main__':unittest.main()
