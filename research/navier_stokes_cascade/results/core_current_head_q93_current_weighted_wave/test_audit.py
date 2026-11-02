import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q93
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q93.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']
class Phase:
    def normal(self,j):return q93.EP(int(j==0))
    def matrix(self,j):return q93.EP(2*int(j==0))
    def moving_normal(self,j):return q93.EP(3*int(j==0))
def make(M=2):
    post=fixture(M);target=post.current.current.state.target
    target.total=lambda j:q93.EP(j+1)
    root=target.root
    return post,q93.CurrentWaves(post,regular_wave=lambda j:3*root.wave(j),phase=Phase(),inverse_normal_square=1)
class Tests(unittest.TestCase):
    def test_current_signed_formula_exact_and_original_root_preserved(self):
        post,w=make();t=w.state.target;root=t.root
        for j in range(3):
            self.assertEqual(w.signed(j,True).record(),q93.convert(t.signed_amplitude(j)).record())
            self.assertEqual(w.signed(j,False).record(),(3*w.signed(j,True)).record())
        self.assertIs(w.cut,root);self.assertIs(w.state.primary.primary.root,root)
        for k in ('H','target','curl','loghalf'):
            self.assertIs(getattr(w.regular,k),getattr(root,k))
        self.assertIsNot(w.regular.wave,root.wave)
    def test_through_norms_distinct_slots_and_pressure(self):
        post,w=make();v=w.through(1)
        root=w.root
        other=q93.CurrentWaves(post,regular_wave=root.wave,phase=Phase(),inverse_normal_square=1).through(1)
        for k in ('T','R','Pi'):self.assertEqual(v[k].record(),(3*other[k]).record())
        self.assertEqual(v['U'].record(),other['U'].record())
        self.assertEqual(v['Pi'].record(),(5*v['T']).record())
        self.assertEqual(v['target_request'],2);self.assertIs(v['wave_state'],w.state)
        self.assertEqual(v['pole'],max(15,12+w.state.target.pole))
        with self.assertRaises(ValueError):w.through(w.state.max_output+1)
    def test_inverse_normal_recurrence_bounds_exact_derivatives(self):
        x=s.symbols('x');_,w=make()
        class Moving(Phase):
            def normal(self,j):return q93.EP(int(j<2))
        w.phase=Moving()
        for j in range(5):
            value=abs(s.diff((1+x)**-2,x,j).subs(x,0))
            expected=w.inverse(j).record()
            # No exponential or S factors occur in this fixture.
            bound=sum(Q(c) for power,terms in expected.items() for c,e in terms)
            self.assertLessEqual(value,bound)
    def test_positive_coarse_reserves_and_label_gaussian_control(self):
        k=s.symbols('k');alpha=2-3*k;B=1-k
        self.assertEqual(s.expand(alpha-k-B),1-3*k)
        for v in (Q(1,100000),Q(1,8),Q(249,1000)):self.assertGreater(1-3*v,0)
        x=s.symbols('x',positive=True)
        # Same-point weights survive raw pullback; distinct Gaussians cannot be identified.
        self.assertEqual(s.limit(s.exp(-x*x)/s.exp(-2*x*x),x,s.oo),s.oo)
        _,w=make();self.assertEqual(w.row(0).record(),q93.EP(5).record())
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
