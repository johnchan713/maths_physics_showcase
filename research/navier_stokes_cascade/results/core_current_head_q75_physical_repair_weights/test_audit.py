import json,unittest
from fractions import Fraction as Q
import sympy as s
import audit,q75

class Tests(unittest.TestCase):
    def test_pullback_and_moving_cutoff(self):
        q,Z=s.symbols('q Z',positive=True);R=q**s.Rational(2,3);D=2/(3*q)
        scale=1+Z; C=3*scale; S=scale**(-s.Rational(3,5))
        y=s.log(R**2/(2*scale)); eta=Z*scale**(-s.Rational(2,5))
        t=s.symbols('t');p=(1+eta**2)*(2+t+t*t);v=eta*t**3
        ref=(2*R*s.diff(S*v.subs(t,y),q),
             (2*C-R**2)*s.diff(S*p.subs(t,y)/R,q)-4*C*D*S*p.subs(t,y)/R)
        result=q75.changes(R,D,S,C,p.subs(t,y),s.diff(p,t).subs(t,y),s.diff(v,t).subs(t,y))
        for a,b in zip(ref,result):
            self.assertEqual(s.simplify(a-b),0)
            self.assertEqual(s.simplify(s.diff(a-b,Z).subs({q:2,Z:1})),0)

    def test_swirl_value_term(self):
        R,D,S,C,P=s.symbols('R D S C P',nonzero=True)
        _,v=q75.changes(R,D,S,C,P,0,0)
        self.assertEqual(s.expand(v-D*S*P*(R**2-6*C)/R),0)

    def test_independent_band_frequency(self):
        # g=1/2: smaller physical bands can fail the direct budget at fixed N.
        for k in (2,4,8):
            eps=Q(1,2**(2*k)); N=4
            lhs=eps**2/N**2;rhs=eps**2*Q(1,2**k)
            self.assertEqual(q75.direct_budget(s.Rational(eps),N,1,1,s.Rational(1,2)),lhs<=rhs)
        self.assertLess(rhs,lhs)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
