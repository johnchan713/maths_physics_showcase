import json,unittest
from fractions import Fraction as Q
import sympy as s
import audit,q67

class Tests(unittest.TestCase):
    def test_forward_memory_normalization(self):
        x,z=s.symbols('x z');g=2+x+3*z+x*z*z
        P=7+z**3+s.integrate(g,(x,1,x))
        for i in range(1,5):
            for j in range(4):
                lhs=s.diff(P,x,i,z,j)/s.factorial(i)/s.factorial(j)
                rhs=s.diff(g,x,i-1,z,j)/s.factorial(i-1)/s.factorial(j)/i
                self.assertEqual(s.expand(lhs-rhs),0)
        self.assertEqual(q67.primitive_bound(0,1,{(0,1):Q(3)},[0,Q(7)],Q(2)),13)
        self.assertEqual(q67.primitive_bound(3,1,{(2,1):Q(7)},[],Q(2)),Q(7,3))

    def test_pressure_memory_and_full_map(self):
        X,E,U,M,Me,S,Se,Pi,Pe,h,z,c,ce=s.symbols('X E U M Me S Se Pi Pe h z c ce')
        D=s.Rational(1,2)-h;A=s.Rational(1,2)+h;d=1-z*z;L=1-2*h*z*z
        W=1-(2*D*z*M+d*Me)/X
        Ns=-W*U+(D*(M-z*Me)+4*h*z*S-d*Se)/X+4*A*z*Pi-d*Pe
        p=X*Ns/(L*E)
        w=X-2*D*z*M-d*Me
        direct=(-w*U+D*(M-z*Me)+4*h*z*S-d*Se+X*(4*A*z*Pi-d*Pe))/(L*E)
        self.assertEqual(s.simplify(p-direct),0)
        change=s.expand(p.subs({Pi:Pi+c,Pe:Pe+ce})-p)
        self.assertEqual(s.simplify(change-X*(4*A*z*c-d*ce)/(L*E)),0)
        k=s.symbols('k',positive=True)
        self.assertEqual(s.simplify(change.subs({z:0,c:0,ce:k})), -k*X/E)

    def test_reciprocal_and_stencil(self):
        x,z=s.symbols('x z');E=2+x+z+x*z;coeff={}
        for total in range(5):
            for i in range(total+1):
                j=total-i
                coeff[i,j]=s.Rational(1,2) if total==0 else -sum(E.coeff(x,b).coeff(z,c)*coeff[i-b,j-c] for b in range(i+1) for c in range(j+1) if b+c>0)/2
        product=s.Poly(s.expand(E*sum(v*x**i*z**j for (i,j),v in coeff.items())),x,z)
        for (i,j),v in product.terms():
            if 0<i+j<=4:self.assertEqual(v,0)
        for m in range(6):
            for l in range(m+1):
                cap=q67.inputs(m,l)
                for i in range(m+1):
                    for j in range(min(l,m-i)+1):
                        self.assertLessEqual(i+1+j,cap['field_total'])
                        self.assertLessEqual(j+1,cap['trace_eta'])

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
