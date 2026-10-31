import json,unittest
import sympy as s
import audit,q59


class Tests(unittest.TestCase):
    def test_squared_solve_with_small_root(self):
        x=s.symbols('x');H=s.Matrix([[1,x],[0,1+x]])
        target=s.Matrix([s.Rational(1,10**12)+x,1+x**2])
        y=H.inv()*target; B=q59.squared(1,[1,1,0,0,0],[1,1,1,0,0])
        for j,upper in enumerate(B):
            self.assertLessEqual(max(abs(a.diff(x,j).subs(x,0)/s.factorial(j)) for a in y),upper)
        self.assertGreater(abs(s.diff(s.sqrt(y[0]),x,2).subs(x,0)),B[2])

    def test_rectangular_ode(self):
        x,v=s.symbols('x v')
        # Scalar solution of t'=(x+v)t, on 0<=x,v<=1, P=1.
        sol=(1+x**2)*s.exp(v*x+v**2/2)
        A=lambda a,b:2 if a+b==0 else int(a+b==1)
        u=q59.path(2,5,5,1,lambda a:2,A)
        self.assertEqual(set(u),{(a,b) for a in range(3) for b in range(6-a)})
        for (a,b),val in u.items():
            self.assertLessEqual(sol.diff(x,a).diff(v,b).subs({x:1,v:1}),val)

    def test_frozen(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
