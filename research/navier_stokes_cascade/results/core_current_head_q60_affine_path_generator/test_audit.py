import json
import unittest

import sympy as s

import audit
from q60 import coefficients, path


class Tests(unittest.TestCase):
    def test_pressure_and_damping_on_rotating_normal(self):
        x, v = s.symbols('x v', real=True)
        n = s.Matrix([x+v, 1, 0])
        m = n.diff(v)
        N = n.dot(n)
        # Nontrivial tangent column, normalized by its value at the base.
        t = s.Matrix([1, -x-v, 0]) / s.sqrt(N) * s.exp(
            -(1+x*x)*v-x*v*v-v**3/3)
        B = -n*m.T - N**2*s.eye(3)
        self.assertEqual(s.simplify(N*t.diff(v)-B*t), s.zeros(3,1))
        self.assertNotEqual((n.T*t.diff(x)).subs({x:0,v:0})[0], 0)
        self.assertNotEqual(((B+n*m.T)*t-N*t.diff(v)).subs({x:0,v:0}), s.zeros(3,1))

        def coeff(p, i, j):
            return p.expand().coeff(x,i).coeff(v,j)
        # The vector 1-norm and matrix induced 1-norm also obey (2).
        def M(i,j):
            return max(sum(abs(coeff(B[r,c],i,j)) for r in range(3)) for c in range(3))
        # Compute the exact Taylor polynomial by total-degree rescaling;
        # no generator recurrence used for the reference values.
        z = s.symbols('z')
        poly = [p.subs({x:z*x,v:z*v}).series(z,0,7).removeO().subs(z,1).expand()
                for p in t]
        def ref(i,j):
            return sum(abs(coeff(p,i,j)) for p in poly)
        bounds = path(3,6,1,lambda a:ref(a,0),lambda i,j:abs(coeff(N,i,j)),M)
        self.assertEqual(len(bounds),22)
        for idx, upper in bounds.items():
            self.assertLessEqual(ref(*idx),upper)

    def test_slow_and_path_coefficient_tables(self):
        x,w = s.symbols('x w')
        base = s.Matrix([1+x, 1-x*x])
        m = s.Matrix([-1+x, x*x])
        K = s.Matrix([[x,1],[-2*x*x,1+x]])
        n = base + w*m
        N = n.dot(n)
        B = -N*K+n*(n.T*K-m.T)-3*N**2*s.eye(2)
        # Entry sums dominate all the real-vector and operator norms.
        def jets(mat):
            return [sum(abs(z.expand().coeff(x,i)) for z in mat) for i in range(5)]
        D,E = coefficients(jets(base),jets(m),jets(K),3)
        for (i,j),upper in D.items():
            self.assertLessEqual(abs(N.expand().coeff(x,i).coeff(w,j)),upper)
        for (i,j),upper in E.items():
            mat = B.applyfunc(lambda z:z.expand().coeff(x,i).coeff(w,j))
            self.assertLessEqual(max(sum(abs(mat[r,c]) for r in range(2))
                                     for c in range(2)),upper)

    def test_same_path_order_denominator_term_is_necessary(self):
        # (1-x)t'=t, t(0)=1. All nonzero rearranged Taylor terms positive.
        tab = path(4,8,1,lambda a:int(a==0),lambda i,j:int((i,j)==(1,0)),
                   lambda i,j:int((i,j)==(0,0)))
        for (a,b),bound in tab.items():
            if b:
                self.assertEqual(bound, s.binomial(a+b-1,a)/s.factorial(b))

    def test_frozen(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))


if __name__ == '__main__':
    unittest.main()
