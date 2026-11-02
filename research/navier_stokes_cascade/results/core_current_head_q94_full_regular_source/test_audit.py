import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q94
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q94.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']
class Phase:
    def normal(self,j):return q94.EP(int(j==0))
    def matrix(self,j):return q94.EP(2*int(j==0))
    def moving_normal(self,j):return q94.EP(3*int(j==0))
class Tests(unittest.TestCase):
    def test_complete_curl_phase_transport_and_zero_output(self):
        R,z=s.symbols('R z',positive=True,real=True);normal=s.Matrix([1,1/R,1])
        C=s.Matrix([R**2,R*z,z**2+R]);D=s.Matrix([R*z+1,R+z,R**2*z])
        def curl(A,m):
            star=s.Matrix([-s.diff(A[1],z),s.diff(A[0],z)-s.diff(A[2],R),s.diff(A[1],R)+A[1]/R])
            return star+s.I*m*normal.cross(A)
        for m,mp in ((1,1),(1,-1),(-2,1),(-1,-2)):
            a,b=curl(C,m),curl(D,mp);div=s.diff(a[0],R)+a[0]/R+s.diff(a[2],z)
            self.assertEqual(s.simplify(div+s.I*m*normal.dot(a)),0)
            self.assertNotEqual(s.simplify(normal.dot(a)),0)
            conn=s.Matrix([-a[1]*b[1],a[1]*b[0],0])/R
            ordinary=a[0]*b.diff(R)+a[2]*b.diff(z)+conn
            direct=ordinary+s.I*mp*a.dot(normal)*b
            reduced=ordinary-s.Rational(mp,m)*div*b
            self.assertEqual((direct-reduced).applyfunc(s.simplify),s.zeros(3,1))
            if mp==-m:self.assertNotEqual(direct.applyfunc(s.simplify),s.zeros(3,1))
    def test_ordered_pair_mass_and_pole_embedding(self):
        # Ordered products include both crosses and the new square.
        old,new=s.symbols('old new');self.assertEqual(s.expand((old+new)**2),old**2+2*old*new+new**2)
        e=s.symbols('e',positive=True);p=3;weight=s.exp(-1/e**2)
        product=weight**2*e**(-2*p)
        target=weight*e**(-2*p)
        self.assertEqual(s.simplify(product/target),weight)
        self.assertEqual(q94.wave_wave(0,1,1,1).record(),q94.EP(6).record())
        with self.assertRaises(ValueError):q94.wave_wave(0,1,1,1,Q(-1))
    def test_typed_full_source_join_and_common_reserve(self):
        post=fixture(2);target=post.current.current.state.target;target.total=lambda j:q94.EP(j+1)
        w=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        source=q94.CurrentSource(w);calls=[];EP=q94.EP
        def vector(name,n):calls.append((name,n));return {k:EP(1) for k in ('x','u','z')}
        post.current.old_full=lambda n:vector('old',n);post.current.fast_unweighted=lambda n:vector('fast',n);post.current.slow=lambda n,b:vector('slow',n)
        post.current.initial.budget._data=lambda n,p:dict(b=EP(1),G=EP(1),c=EP(1),graph=EP(2),r=EP(3))
        v=source.coefficient(0,{k:EP(1) for k in ('Eik','bR','Bz','psi')})
        self.assertEqual(v['F'].record(),(v['linear']+v['mean_wave']+v['wave_wave']).record())
        self.assertEqual(calls,[('old',2),('fast',2),('slow',2)])
        self.assertEqual(v['common_pole'],2*w.through(2)['pole']);self.assertEqual(v['mean_factory_reserve'],2)
        self.assertFalse(v['actual_native_contracts_certified']);self.assertTrue(v['wave_wave'].record())
        with self.assertRaises(ValueError):source.coefficient(1,{k:EP(1) for k in ('Eik','bR','Bz','psi')})
        with self.assertRaises(ValueError):source.coefficient(0,{})
        with self.assertRaises(TypeError):q94.CurrentSource(lambda n:1)
    def test_input_ratio_not_output_frequency(self):
        # m=1,m'=-1 has rho=1 and is well-defined despite zero output.
        self.assertEqual(q94.wave_wave(2,3,5,7,1).record(),q94.EP(2*(3*(3+1)+5)*49).record())
        self.assertEqual(q94.record()['source_order'],'1-kappa')
        self.assertFalse(q94.record()['output_harmonic_division'])
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
