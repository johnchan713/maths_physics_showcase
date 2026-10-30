from fractions import Fraction as Q
import unittest,json
import sympy as s
import completion as c
import audit

class CompletionTests(unittest.TestCase):
    def test_threshold_slow_and_length_gates(self):
        for h,r,L in ((Q(1,200),Q(1,6057072),Q(100)),(Q(1,10**20),Q(1,100),Q(10**120))):
            n=c.threshold(h,r,L)
            self.assertGreaterEqual(n,16)
            self.assertLessEqual(360*(4/h)**5/n,1)
            self.assertGreaterEqual(2*r*n*n,L)
            self.assertLessEqual(12*r,1)

    def test_rounding_floor_for_nonzero_signed_carriers(self):
        for k in (1,3,1000):
            for j in (-7,-1,1,9):
                for R in (Q(1,2),Q(2),Q(5)):
                    self.assertGreaterEqual(abs(Q(j,k)/R),Q(1,5*k))

    def test_zero_carrier_would_destroy_floor(self):
        self.assertLess(Q(0),Q(1,5000))

    def test_finite_head_budget_and_union(self):
        h=c.head_recipe(100,5,12,3,2,7,4)
        self.assertEqual(h['normal_floor'],Q(1,500))
        self.assertEqual(h['log_propagator_constant'],12*(6+1000+7+4))
        u=c.combine(h,10,4,Q(1,4))
        self.assertEqual(u['log_T'],h['log_propagator_constant'])
        self.assertEqual(u['normal_reciprocal'],250000)
        self.assertEqual(u['damping_floor'],Q(1,250000))

    def test_empty_head(self):
        self.assertEqual(c.combine(None,10,4,Q(1,4)),{'log_T':10,'normal_reciprocal':4,'damping_floor':Q(1,4)})

    def test_constraint_with_moving_normal(self):
        v=s.symbols('v',real=True);n=s.Matrix([v,1,2]);np=n.diff(v)
        K=s.Matrix([[0,-2,0],[3,0,0],[4,0,0]]);d=2
        A=-K+n*(n.T*K-np.T)/(n.dot(n))-d*s.eye(3)
        self.assertEqual((np.T+n.T*A+d*n.T).applyfunc(s.simplify),s.zeros(1,3))

    def test_noncommuting_scalar_damping_identity(self):
        v=s.symbols('v',real=True);m=s.symbols('m',integer=True)
        V=s.Matrix([[1,v],[v,1+v*v]])
        A=V.diff(v)*V.inv();D=v+2;I=v*v/2+2*v
        Vm=s.exp(-(m*m-1)*I)*V
        self.assertEqual((Vm.diff(v)-(A-(m*m-1)*D*s.eye(2))*Vm).applyfunc(s.simplify),s.zeros(2))

    def test_upper_bound_does_not_give_component_floor(self):
        v=s.symbols('v',real=True);V=s.Matrix([[s.cos(v),-s.sin(v)],[s.sin(v),s.cos(v)]])
        self.assertEqual((V.T*V).applyfunc(s.simplify),s.eye(2))
        self.assertEqual(V[0,0].subs(v,s.pi/2),0)
        self.assertEqual(V.det().simplify(),1)

    def test_guards(self):
        with self.assertRaises(ValueError):c.threshold(Q(1,200),1,1)
        with self.assertRaises(TypeError):c.threshold(.005,Q(1,100),1)
        with self.assertRaises(ValueError):c.head_recipe(0,1,1,1,1,1,1)
        with self.assertRaises(ValueError):c.combine(None,-1,1,1)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['root_or_cone_lower_bound_proved'])

if __name__=='__main__':unittest.main()
