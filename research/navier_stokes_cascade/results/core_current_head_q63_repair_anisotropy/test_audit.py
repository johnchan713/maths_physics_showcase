import json, unittest
from math import comb,factorial
import sympy as s
import audit,q63


class Tests(unittest.TestCase):
    def test_exact_flat_polynomial_majorant(self):
        u=s.symbols('u');p=s.Integer(1)
        for j in range(14):
            terms=s.Poly(p,u).terms()
            old=sum(abs(c)*factorial((int(m[0])+1)//2) for m,c in terms)
            self.assertLessEqual(old,20**j*factorial(j)**3)
            p=s.expand(-u*u*s.diff(p,u)+2*u**3*p)
        b=[81]
        for j in range(1,15):
            b.append(162*sum(20**k*b[j-k] for k in range(1,j+1)))
            self.assertEqual(b[j],81*162*20**j*163**(j-1))
            self.assertLessEqual(b[j],81*3260**j)

    def test_same_root_and_separation(self):
        y,t=s.symbols('y t')
        c=(s.sqrt(1+4*t)-1)/2
        bounds=q63.repair(1,1,[0,1,0,0,0,0,0])
        e=1+t; beta=y**7
        for k in range(7):
            C=bounds[k]+(bounds[k-1] if k else 0)
            for r in range(8):
                exact=s.diff(e*c*beta,t,k,y,r).subs({t:0,y:1})
                self.assertLessEqual(abs(exact),factorial(k)*C*factorial(7)/factorial(7-r))

    def test_chain_rule_eta_count(self):
        R,Z,T=s.symbols('R Z T')
        X=R**2/(2*(1+T));eta=Z/(1+T)**s.Rational(2,5)
        f=s.Function('f')
        expr=s.diff(f(X,eta),R,5,Z,1,T,1)
        # Every derivative of the outer function has eta order at most two.
        derivatives=expr.atoms(s.Derivative)
        self.assertTrue(derivatives)
        for derivative in derivatives:
            counts=dict(derivative.variable_count)
            eta_arg=derivative.expr.args[1]
            self.assertLessEqual(counts.get(eta_arg,0),2)
        self.assertEqual(max(dict(d.variable_count).get(d.expr.args[1],0)
                             for d in derivatives),2)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))


if __name__=='__main__':unittest.main()
