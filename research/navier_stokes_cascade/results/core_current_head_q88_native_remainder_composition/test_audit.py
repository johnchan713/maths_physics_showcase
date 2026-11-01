import unittest,json,importlib.util
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
import audit,q88
p=Path(__file__).resolve().parent.parent
sp=importlib.util.spec_from_file_location('q88_old',p/'core_current_head_q84_nonlinear_radial_envelope/q84.py');old=importlib.util.module_from_spec(sp);sp.loader.exec_module(old)
sp=importlib.util.spec_from_file_location('q88_fixture',p/'core_current_head_q87_final_angular_adapter/test_audit.py')
# Use that fixture's source with this module instance's exact class tree.
fixture_source=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text()
start=fixture_source.index('def fixture(');end=fixture_source.index('class Tests',start)
ns={'q87':q88.angular,'Q':Q};exec(fixture_source[start:end],ns);fixture=ns['fixture']
class Tests(unittest.TestCase):
    def test_same_polynomial_as_frozen_rational_control(self):
        for m in (0,1,3):
            data={k:Q(i+1,3) for i,k in enumerate(('b0','V','G','x','u','z','W','Z','e','r','g','c','T'))}
            self.assertEqual(q88.polynomial(m,data).record(),q88.EP(old.envelope(m,data)['projected_remainder']).record())
    def test_temporal_quadratic_norm_and_shell_count(self):
        a,b=s.symbols('a b',integer=True)
        product=(a+(1+s.sqrt(2))*b)*(a+(1-s.sqrt(2))*b)
        self.assertEqual(s.expand(product-(a*a+2*a*b-b*b)),0)
        for K in (1,2,4,8):
            modes=[(i,j) for i in range(-2*K,2*K+1) for j in range(-2*K,2*K+1) if K<=abs(i)+abs(j)<2*K]
            self.assertLessEqual(len(modes),25*K*K)
            for i,j in modes:
                norm=i*i+2*i*j-j*j
                self.assertNotEqual(norm,0)
                self.assertLessEqual(abs(i)+3*abs(j),6*K)
        self.assertEqual(25*6*16*2*2,9600)
    def test_native_geometry(self):
        U=s.symbols('U',positive=True)
        for n in range(4):
            g=q88.geometry(n,1,4)
            norm=lambda f:sum(abs(s.diff(f,U,j).subs(U,1))/s.factorial(j) for j in range(n+1))
            self.assertLessEqual(norm(U**-2),g['r']);self.assertLessEqual(norm(1/(2*U)),g['g'])
        self.assertEqual(q88.geometry(0,Q(1,2),4)['Umin'],Q(1,8))
    def test_typed_native_composition_manufactured_leaves(self):
        post=fixture(7);out=q88.NativeRemainder(post);f=out.means;EP=q88.EP
        vector=lambda n:{k:EP(1) for k in ('x','u','z')}
        post.current.old_full=vector;post.current.fast_unweighted=vector;post.current.slow=lambda n,w:vector(n)
        budget=post.current.initial.budget
        budget.fixed=lambda n,p:{k:EP(1) for k in ('b','V','G','sigma','c','rho','Kp','Kmom')}
        budget._data=lambda n,p:{k:EP(1) for k in ('b','V','G','graph','r')}
        budget.compact=lambda n,p:dict(KT=EP(1),KA=EP(1),I=EP(1))
        f.state.target.source.haar.rows=lambda n:dict(Etheta=EP(1))
        f.state.coefficients=lambda n:{k:EP(1) for k in ('A','B','L','S')}
        f.state.covariance=lambda n:dict(complete=EP(7),pole=20)
        v=out.inputs(0)
        self.assertEqual(v['chart'],'factorial-through-native-U');self.assertEqual(v['factory_reserve'],7)
        self.assertEqual(v['data']['T'].record(),EP({1:9600*24}).record())
        self.assertTrue(q88.polynomial(0,v['data']).record())
        with self.assertRaises(ValueError):out.inputs(1)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
