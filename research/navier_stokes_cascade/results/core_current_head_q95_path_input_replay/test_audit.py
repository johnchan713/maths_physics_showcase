import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q95
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q95.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']
class Phase:
    def normal(self,j):return q95.EP(int(j==0))
    def matrix(self,j):return q95.EP(2*int(j==0))
    def moving_normal(self,j):return q95.EP(3*int(j==0))
def constant(v):
    rec=v.record();assert set(rec)<={'0'}
    return sum(Q(c) for terms in rec.values() for c,e in terms if e=='0')
class Tests(unittest.TestCase):
    def test_mixed_Gaussian_quotient_majorant(self):
        x,y=s.symbols('x y');P=s.exp(-x*x-y*y);u=1+x+y+x*y;f=P*u
        def maximum(z,j):return max(abs(s.diff(z,x,i,y,j-i).subs({x:0,y:0})) for i in range(j+1))
        raw=[Q(maximum(f,j)) for j in range(4)]
        log=[Q(12345)]+[Q(maximum(-x*x-y*y,j)) for j in range(1,4)]
        out=q95.quotient(raw,log,[0]*4)
        for j in range(4):self.assertLessEqual(maximum(u,j),constant(out['jets'][j]))
        other=q95.quotient(raw,[0]+log[1:],[0]*4)
        self.assertEqual([z.record() for z in out['jets']],[z.record() for z in other['jets']])
    def test_log_poles_are_accumulated(self):
        out=q95.quotient([1]*4,[0,1,1,1],[0,2,1,0])
        self.assertEqual(out['Bell_poles'],[0,2,4,6]);self.assertEqual(out['extra_pole'],6)
        self.assertEqual(q95.quotient([1]*4,[0,1,1,1],[0]*4)['extra_pole'],0)
        with self.assertRaises(ValueError):q95.quotient([1],[0],[Q(1,2)])
        with self.assertRaises(ValueError):q95.quotient([1],[0,1],[0])
    def test_typed_source_path_normal_and_six_variable_count(self):
        post=fixture(3);target=post.current.current.state.target;target.total=lambda j:q95.EP(j+1);EP=q95.EP
        w=q95.source.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        source=q95.source.CurrentSource(w);out=q95.PathInputs(source)
        vector=lambda n:{k:EP(1) for k in ('x','u','z')}
        post.current.old_full=vector;post.current.fast_unweighted=vector;post.current.slow=lambda n,b:vector(n)
        post.current.initial.budget._data=lambda n,p:dict(b=EP(1),G=EP(1),c=EP(1),graph=EP(2),r=EP(3))
        table=[{k:EP(1) for k in ('Eik','bR','Bz','psi')} for j in range(2)]
        v=out.build(1,table,[EP(0),EP(1)],[0,2]);G=v['normalized']
        self.assertEqual(v['source_pole'],max(source.coefficient(j,table[j])['common_pole'] for j in range(2)))
        self.assertEqual(v['common_pole'],v['source_pole']+2);self.assertEqual(v['source_factory_reserve'],3)
        self.assertEqual(v['raw_path'][1].record(),(2*source.coefficient(1,table[1])['F']).record())
        self.assertEqual(v['six_variable_through'].record(),(G[0]+6*G[1]).record())
        self.assertNotEqual(v['six_variable_through'].record(),(G[0]+5*G[1]).record())
        self.assertTrue(v['source_normal_pressure'][0].record());self.assertFalse(v['endpoint_flatness_from_norm_alone'])
        with self.assertRaises(ValueError):out.build(1,table[:1],[0,1])
        with self.assertRaises(TypeError):q95.PathInputs(lambda n:1)
    def test_norm_bound_does_not_give_entrance_flatness_or_zero_pressure(self):
        v=s.symbols('v');P=s.exp(-(v-1)**2);f=P
        self.assertEqual(s.simplify(f/P),1);self.assertNotEqual(f.subs(v,0),0)
        n=s.Matrix([1,0,0]);self.assertEqual((s.eye(3)-n*n.T)*n,s.zeros(3,1))
        pi=s.I*f/2;self.assertEqual(s.simplify(pi/P),s.I/2)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
