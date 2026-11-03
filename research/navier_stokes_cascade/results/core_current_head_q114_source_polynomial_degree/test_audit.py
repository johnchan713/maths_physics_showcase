import unittest,json,math
from fractions import Fraction as Q
from pathlib import Path
import audit,q114
E=q114.radial.EP
text=(audit.HERE.parent/'core_current_head_q110_radial_tail_request/test_audit.py').read_text()
ns=dict(q110=q114.radial,audit=audit,Path=Path,Q=Q)
exec(text[text.index('def setup('):text.index('\nclass Tests(')],ns);setup=ns['setup']

class Tests(unittest.TestCase):
    def test_zero_positive_algebra_and_fixed_exponent(self):
        z=q114.as_growth(0);a=q114.Growth(3,4);b=q114.Growth(2,9)
        self.assertEqual(z+a,a);self.assertEqual((z*a).degree,None)
        self.assertEqual(a+b,q114.Growth(3,10));self.assertEqual(a*b,q114.Growth(5,13))
        p=E({2:3})
        self.assertEqual(q114.compress(p).degree,2)
        # Obtain the original fixed-coefficient ring without a second module identity.
        cls=type(next(iter(E(1).c.values())))
        self.assertEqual(q114.compress(E({2:cls([(2,100)])})),q114.Growth(2,102))

    def test_degree_matches_independent_current_source_expansion(self):
        for n in range(4):
            for A in (0,2):
                mean={k:E({j:2}) for k,j in zip(('x','u','z'),(1,3,2))}
                fixed={k:E({j:3}) for k,j in zip(('b','V','G','c','graph'),(2,1,4,2,9))}
                W=E({4:5});absorb=E({A:2})
                source=q114.radial.radial.radial_envelope(n,mean,W,Q(1,2),fixed,absorb)['reduced_source']
                bound=q114.reduced(n,{k:q114.compress(v) for k,v in mean.items()},q114.compress(W),Q(1,2),{k:q114.compress(v) for k,v in fixed.items()},q114.compress(absorb))
                x,u,z=1+A,3+A,2+A
                expected=max(2+x,2+x,2*x,4,2+z,4+x,x+z,1+u,2*u)
                self.assertEqual(max(source.c),expected);self.assertEqual(bound.degree,expected)
                for S in (0,1,3,10):
                    terms=[math.log(float(c))+float(e)+(k*math.log(S) if k else 0) for k,v in source.c.items() if S or k==0 for e,c in v.terms.items()]
                    if not terms:continue # source vanishes at S=0
                    m=max(terms);actual=m+math.log(sum(math.exp(t-m) for t in terms))
                    self.assertLessEqual(actual,float(bound.height)+bound.degree*math.log1p(S))
                alt=dict(fixed,graph=E({100:999}))
                other=q114.radial.radial.radial_envelope(n,mean,W,Q(1,2),alt,absorb)['reduced_source']
                self.assertEqual(source.record(),other.record())

    def test_all_components_and_covariance_survive(self):
        zero={k:0 for k in ('x','u','z')};fixed={k:1 for k in ('b','V','G','c','graph')}
        self.assertIsNone(q114.reduced(0,zero,0,1,fixed,1).degree)
        for key in zero:
            self.assertEqual(q114.reduced(0,dict(zero,**{key:q114.Growth(5)}),0,1,fixed,1).degree,5 if key=='z' else 10)
        self.assertEqual(q114.reduced(0,zero,q114.Growth(7),1,fixed,1).degree,7)

    def test_current_source_route_and_archive(self):
        obj,b=setup();scope=q114.scope.CompactShapeScope(q114.scope.profiles.KernelProfiles(q114.scope.profiles.charts.ChartBridge(obj)))
        proof=q114.SourceGrowth(scope);rad=obj.radial
        rad.fields=lambda n,b:dict(fields=dict(x=E(1),u=E(2),z=E(3)),pole=b['covariance_pole'])
        rad.current.initial.budget._data=lambda n,p:{k:E(1) for k in ('b','V','G','c','graph')}
        # Actual graph routing, with manufactured norm leaves only.
        out=proof.source(0,b)
        self.assertIs(out['source']['source_state'],obj.state)
        self.assertEqual(out['growth'].degree,0)
        self.assertTrue(out['native_leaf_inequalities_still_required']);self.assertFalse(out['actual_A0_gain_certified'])
        with self.assertRaises(ValueError):proof.source(100005,b)
        obj.state=object()
        with self.assertRaises(ValueError):proof.source(0,b)

    def test_invalid_inputs(self):
        for args in ((-1,0),(1,-1)):
            with self.assertRaises(ValueError):q114.Growth(*args)
        with self.assertRaises(ValueError):q114.as_growth(-1)
        with self.assertRaises(TypeError):q114.SourceGrowth(None)

    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
