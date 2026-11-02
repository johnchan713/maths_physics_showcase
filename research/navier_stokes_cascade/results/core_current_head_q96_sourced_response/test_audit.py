import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q96
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q96.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']
class Phase:
    def normal(self,j):return q96.EP(int(j==0))
    def matrix(self,j):return q96.EP(2*int(j==0))
    def moving_normal(self,j):return q96.EP(3*int(j==0))
class Tests(unittest.TestCase):
    def test_normalized_ambient_and_full_forced_pressure_exact(self):
        v,x=s.symbols('v x',real=True);n=s.Matrix([1,v+x,0]);np=n.diff(v);N=n.dot(n)
        K=s.Matrix([[0,2,0],[3,0,1],[1,0,0]]);P=s.exp(-v*v);rate=-2*v
        H=n*n.T/N;A0=-K+n*(n.T*K-np.T)/N
        C=A0-(N+rate)*s.eye(3)+(N+rate)*H
        self.assertEqual((np.T+n.T*C).applyfunc(s.simplify),s.zeros(1,3))
        self.assertNotEqual((np.T+n.T*(C+rate*s.eye(3))).applyfunc(s.simplify),s.zeros(1,3))
        a=v*v*(1+x);u=s.Matrix([-(v+x)*a,a,x*v*v]);g=C*u-u.diff(v)+n*(1+x)
        self.assertEqual(s.simplify(n.dot(u)),0);self.assertNotEqual(s.simplify(n.dot(u.diff(x))),0)
        pi=P*s.I*((K.T*n-np).dot(u)+n.dot(g))/(2*N)
        t=P*u;f=P*g
        residual=t.diff(v)+K*t+N*t+2*s.I*n*pi+f
        self.assertEqual(residual.applyfunc(s.simplify),s.zeros(3,1))
        self.assertEqual(t.subs(v,0),s.zeros(3,1))
    def test_normalized_generator_and_triangular_budget(self):
        v=q96.response([1,0],[1,0],[0,0],[0,0],[1,0],1,0,1)
        self.assertEqual(v['generator'][0].record(),q96.EP(10).record())
        self.assertEqual(v['velocity'][0].record(),q96.EP(2).record())
        self.assertEqual(v['velocity'][1].record(),q96.EP(21).record())
        self.assertEqual(v['pressure'][0].record(),q96.EP(1).record())
        with self.assertRaises(ValueError):q96.response([1],[1],[0],[0],[0],1,1,1)
    def test_typed_cut_response_evaluation_and_cutoff_curl(self):
        post=fixture(3);target=post.current.current.state.target;target.total=lambda j:q96.EP(j+1);EP=q96.EP
        q94=q96.inputs.source
        w=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        out=q96.SourcedResponse(q96.inputs.PathInputs(q94.CurrentSource(w)))
        vector=lambda n:{k:EP(1) for k in ('x','u','z')}
        post.current.old_full=vector;post.current.fast_unweighted=vector;post.current.slow=lambda n,b:vector(n)
        budget=post.current.initial.budget
        budget._data=lambda n,p:dict(b=EP(1),G=EP(1),c=EP(1),graph=EP(2),r=EP(3))
        budget.compact=lambda n,p:dict(KT=EP(1),KA=EP(1),I=EP(1))
        table=[{k:EP(1) for k in ('Eik','bR','Bz','psi')} for j in range(2)]
        v=out.build(0,table,[0,1],[1,0],[1,0],tangent=EP(1),length=EP(1),log_poles=[0,2])
        changed=out.build(0,table,[0,1],[1,0],[1,1],tangent=EP(1),length=EP(1),log_poles=[0,2])
        self.assertEqual(v['evaluation_factor'].record(),EP({0:1,1:376}).record())
        self.assertEqual(v['principal'][0].record(),changed['principal'][0].record())
        self.assertNotEqual(v['curl'][0].record(),changed['curl'][0].record())
        self.assertEqual(v['complete'][0].record(),(v['principal'][0]+v['curl'][0]).record())
        self.assertEqual(v['source_ceiling'],1);self.assertEqual(v['source_factory_reserve'],3)
        base=out.inputs.build(1,table,[0,1],[0,0])['source_pole'];self.assertEqual(v['pole'],base+4)
        self.assertFalse(v['actual_contracts_certified'])
        with self.assertRaises(TypeError):q96.SourcedResponse(lambda n:1)
    def test_normal_only_source_retains_pressure_and_exit_need_not_vanish(self):
        v=s.symbols('v');n=s.Matrix([1,0,0]);f=n*v**2*(1-v)**2
        self.assertEqual((s.eye(3)-n*n.T)*f,s.zeros(3,1));self.assertNotEqual(s.I*n.dot(f)/2,0)
        self.assertNotEqual(s.integrate(v**2*(1-v)**2,(v,0,1)),0)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
