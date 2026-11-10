import json
import tempfile
import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q140 as q
import audit

class Tests(unittest.TestCase):
    def fixture(self):
        return q.four_to_six(4,2,3,4,2,1,1,2,8,1,q.IDENTITY)

    def test_four_to_six_and_q139_binding(self):
        r=self.fixture()
        self.assertEqual(set(r['components']),q.KINDS)
        self.assertGreaterEqual(r['H'],sum(r['components'].values()))
        disk=q.previous.phase_disk(q.previous.scalar.parameter_disk(2,8,q.previous.scalar.IDENTITY),q.previous.IDENTITY)
        joined=q.bind(r,disk,1,8,q.IDENTITY)
        self.assertEqual(joined['primitive']['scalar_jet_prefactor'],r['H'])
        self.assertEqual(joined['pullback']['N'],8)
        self.assertFalse(joined['actual_regional_four_field_mixed_bound_certified'])

    def test_finite_composition_counting_inequality(self):
        # n! prod(i_l! k_l!) <= (i+k)! <= 2^(i+k)i!k!.
        from itertools import product
        # A small exhaustive grid checks the combinatorics without turning
        # the provenance suite into an exponential stress test.
        for i in range(4):
            for k in range(4):
                m=i+k
                for n in range(1,m+1):
                    for ys in product(range(i+1),repeat=n):
                        if sum(ys)!=i: continue
                        for es in product(range(k+1),repeat=n):
                            if sum(es)!=k or any(ys[j]+es[j]==0 for j in range(n)): continue
                            lhs=factorial(n)
                            for a,b in zip(ys,es): lhs*=factorial(a)*factorial(b)
                            self.assertLessEqual(lhs,2**m*factorial(i)*factorial(k))

    def test_flat_cutoff_and_augmented_variance_controls(self):
        import sympy as s
        a,b,z=s.symbols('a b z', positive=True)
        v0=a+b*b/a; v=v0+a*z*z
        self.assertEqual(s.simplify(v-v0-a*z*z),0)
        r=self.fixture()
        self.assertGreaterEqual(r['components']['v'],r['unaugmented_variance_cap'])
        self.assertFalse(r['flat_cutoff_slow_analytic_extension_claimed'])

    def test_inherited_q72_cap_obstruction(self):
        # A modest order already separates the inherited cap from a very large
        # manufactured fixed Gevrey-3 budget.  The README proves eventual
        # separation for every fixed pair of constants.
        lower=q.q72_log2_lower(Q(1,100),Q(1,100),6)
        upper=q.gevrey3_log2_budget(10**30,10**30,6)
        self.assertGreater(lower,upper)
        self.assertGreater(q.q72_log2_lower(Q(1,100),Q(1,100),7),lower)

    def test_bounds_binary_and_guards(self):
        r=self.fixture()
        for kind in q.KINDS:
            for i in range(5):
                for k in range(5):
                    self.assertLessEqual(q.bound(r,kind,i,k)['upper'],2**q.bound(r,kind,i,k,True)['upper'])
            self.assertIsInstance(q.bound(r,kind,100005,100005,True)['upper'],int)
        with self.assertRaises(ValueError): q.four_to_six(4,2,3,4,2,1,1,2,8,1,'different-cutoff')
        with self.assertRaises(ValueError): q.bound(r,'s',0,0)
        with self.assertRaises(ValueError): q.q72_log2_lower(1,1,0)

    def test_00_record(self):
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/'replay.json'; frozen=audit.HERE/'evidence.json'
            out.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+'\n')
            self.assertEqual(out.read_bytes(),frozen.read_bytes())

if __name__ == '__main__': unittest.main()
