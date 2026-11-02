import unittest,json
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
import audit,q91
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q91.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']
class Tests(unittest.TestCase):
    def test_exact_cartesian_transport_in_cylindrical_basis(self):
        R,z,th=s.symbols('R z theta',positive=True,real=True)
        er=s.Matrix([s.cos(th),s.sin(th),0]);et=s.Matrix([-s.sin(th),s.cos(th),0]);ez=s.Matrix([0,0,1]);E=s.Matrix.hstack(er,et,ez)
        m=s.Matrix([R,z,R*z]);a=s.Matrix([R**2,z**2,R*z+1]);phase=th+R+z;wave=a*s.exp(s.I*phase)
        def adv(u,v):
            cart=E*v
            return u[0]*cart.diff(R)+u[1]/R*cart.diff(th)+u[2]*cart.diff(z)
        actual=E.T*(adv(m,wave)+adv(wave,m))*s.exp(-s.I*phase)
        n=s.Matrix([1,1/R,1]);conn=s.Matrix([-2*m[1]*a[1],m[1]*a[0]+m[0]*a[1],0])/R
        predicted=m[0]*a.diff(R)+m[2]*a.diff(z)+a[0]*m.diff(R)+a[2]*m.diff(z)+s.I*m.dot(n)*a+conn
        self.assertEqual((actual-predicted).applyfunc(s.simplify),s.zeros(3,1))
        self.assertNotEqual(conn,s.zeros(3,1));self.assertNotEqual(s.simplify(m.dot(n)),0)
    def test_weight_is_supplied_by_wave_not_mean(self):
        x,y=s.symbols('x y',positive=True);weight=s.exp(-1/x**2)
        self.assertEqual(s.simplify((weight*s.sin(y))/weight),s.sin(y))
        self.assertEqual(s.limit(1/weight,x,0,dir='+'),s.oo)
        # Product-rule control with raw jets; no derivative of the weight is deleted.
        mean=1+x;wave=weight*s.sin(y)
        self.assertEqual(s.simplify(s.diff(mean*wave,x)-s.diff(mean,x)*wave-mean*s.diff(wave,x)),0)
    def test_factorial_derivative_shift(self):
        x,y=s.symbols('x y');f=(1+x+y)**5
        def norm(g,n):return sum(abs(s.diff(g,x,i,y,j-i).subs({x:0,y:0}))/s.factorial(i)/s.factorial(j-i) for j in range(n+1) for i in range(j+1))
        for n in range(4):self.assertLessEqual(norm(s.diff(f,x),n),(n+1)*norm(f,n+1))
    def test_typed_mean_composition_and_state_guard(self):
        post=fixture(1);out=q91.MeanWave(post);calls=[];EP=q91.EP
        def vec(name,n):calls.append((name,n));return {k:EP(i+1) for i,k in enumerate(('x','u','z'))}
        post.current.old_full=lambda n:vec('old',n);post.current.fast_unweighted=lambda n:vec('fast',n);post.current.slow=lambda n,w:vec('slow',n)
        post.current.initial.budget._data=lambda n,p:dict(graph=EP(2),r=EP(3))
        v=out.coefficient(0,wave=EP(5),normal=EP(7),wave_state=out.means.state)
        self.assertEqual(calls,[('old',1),('fast',1),('slow',1)])
        # M=243*3*(1+2+3), bracket=2*1*3+2*7+2*3=26.
        self.assertEqual(v['coefficient'].record(),EP(3*26*243*18*5).record())
        self.assertEqual(v['field_and_wave_order'],1);self.assertEqual(v['mean_weight_required'],0)
        with self.assertRaises(ValueError):out.coefficient(0,wave=EP(5),normal=EP(7),wave_state=object())
        with self.assertRaises(ValueError):out.coefficient(1,wave=EP(5),normal=EP(7),wave_state=out.means.state)
        with self.assertRaises(TypeError):q91.MeanWave(lambda n:1)
        with self.assertRaises(ValueError):q91.envelope(0,1,1,1,1,1,harmonic=0)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
