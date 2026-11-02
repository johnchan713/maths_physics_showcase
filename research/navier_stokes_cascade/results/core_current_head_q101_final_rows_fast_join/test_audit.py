import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q101
q100=q101.changes;q99=q100.scalars;q98=q100.pressure
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q101.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

class Tests(unittest.TestCase):
    def test_weight_absorption_and_pole_embedding(self):
        for w in (Q(1,10),Q(1),Q(4)):
            for degree in (1,2):self.assertGreaterEqual(w*(5*3)**degree,(w*3)**degree)
        self.assertLess(4*3**2,(4*3)**2) # naive same-pole squaring fails
        self.assertEqual(q101.embedding(Q(3),2,5),27)
        self.assertEqual(q101.embedding(Q(1,2),2,5),1)
        with self.assertRaises(ValueError):q101.embedding(3,5,2)
    def test_exact_final_cancellation_keeps_postsigned_tail(self):
        y=s.symbols('y',real=True);rf=s.sin(y);rs=2*s.cos(2*y)
        beforet=3+s.cos(y);beforez=5+2*s.cos(y)+s.sin(2*y)
        fastt=-s.cos(y);fastz=-s.sin(2*y)
        finalt=beforet+fastt+rf+rs;finalz=beforez+fastz+rf+rs
        Qcirc=lambda f:s.expand(f-s.integrate(f,(y,0,2*s.pi))/(2*s.pi))
        self.assertEqual(s.simplify(Qcirc(finalt)-(rf+rs)),0)
        self.assertEqual(s.simplify(Qcirc(finalz)-(2*s.cos(y)+rf+rs)),0)
        self.assertNotEqual(s.simplify(Qcirc(finalz)-(s.sin(3*y)+rf+rs)),0)
    def test_positive_total_join(self):
        E=q101.EP
        v=q101.join(dict(Etheta=E(2),Ez=E(3)),dict(Etheta=E(5),Ez=E(7)),11,dict(Qtheta=E(13),Qz=E(17)),2,3)
        self.assertEqual(v['incoming_theta'].record(),E(14).record())
        self.assertEqual(v['incoming_z'].record(),E(64).record())
        self.assertEqual(v['Etheta'].record(),E(53).record())
        self.assertEqual(v['Ez'].record(),E(115).record())
    def test_typed_join_domains_and_complete_fast_inverse(self):
        q97=q98.covariance;q94=q97.response.inputs.source;post=fixture(8);E=q101.EP
        class Phase:
            def normal(self,j):return E(int(j==0))
            def matrix(self,j):return E(2*int(j==0))
            def moving_normal(self,j):return E(3*int(j==0))
        target=post.current.current.state.target
        wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        cov=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
        oldpole=post.current.fast.rows.pole;pole=oldpole+1
        pm=q98.ShiftedPressure(cov,c0=2,c1=4,shape=[1]*8,time_inverse=lambda n,p:1,maximum=6,pole=pole)
        change=q100.SourcedRowChanges(q99.ChangedScalars(pm))
        adapter=q101.FinalRows(change,maximum=5,pole=pole,c0=2,c1=4,shape=iter([1]*8),time_inverse=lambda n,p:1)
        calls=[];vector=lambda a:{k:E(a) for k in ('x','u','z')}
        post.current.old_full=lambda n:vector(1)
        post.current.slow=lambda n,w:vector(2 if w else 3)
        post.current.initial.budget._data=lambda n,p:{k:E(1) for k in ('b','V','G','c','graph','r')}
        post.current.fast.rows.rows=lambda n:(calls.append(('postsigned',n)) or dict(Etheta=E(3),Ez=E(5)))
        post.current.fast.compact=lambda n,p:(calls.append(('oldA1',n,p)) or dict(order=n,pole=p,KT=E(2),KA=E(3),I=E(7)))
        bundle=dict(source_state=cov.state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=2,covariance_pole=pole,change=[E(1)]*3)
        # Real new T0 recipes and current pressure are used at output zero;
        # old coefficient leaves are manufactured, not native certificates.
        out=adapter.bound(0,bundle)
        self.assertEqual(out['pole'],pole);self.assertFalse(out['means_updated'])
        self.assertIn(('oldA1',0,oldpole),calls)
        self.assertEqual(out['post_signed_A1'].record(),E(3*243*5).record())
        with self.assertRaises(ValueError):adapter.bound(0,dict(bundle,covariance_pole=pole+1))
        with self.assertRaises(ValueError):adapter.bound(0,dict(bundle,means_updated=True))
        # Isolate the inherited inverse interface from enormous nested leaves.
        seen=[]
        adapter.bound=lambda n,b:(seen.append(n) or dict(Etheta=E(n+1),Ez=E(n+2)))
        adapter.axial=lambda n,p:dict(order=n,pole=p,KT=E(3),KA=E(5),I=E(7))
        fast=adapter.fast(0,bundle)
        self.assertEqual(seen,[4,5]);self.assertEqual((5*fast['z']).record(),(36*fast['u']).record())
        self.assertEqual(fast['u'].record(),(E({1:322560})).record())
        self.assertFalse(fast['fresh_slow_or_postfast_scalars_supplied'])
        with self.assertRaises(ValueError):adapter.fast(1,bundle)
        with self.assertRaises(TypeError):q101.FinalRows(None,maximum=0,pole=1,c0=2,c1=4,shape=[1]*3,time_inverse=lambda n,p:1)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
