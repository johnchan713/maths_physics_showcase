import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q108
q105=q108.projected;q104=q105.radial;q103=q104.post;q102=q103.slow;q101=q102.rows
q100=q101.changes;q99=q100.scalars;q98=q100.pressure
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text()
ns={'q87':type('Module',(),{'full':q108.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

def setup():
    post=fixture(13);E=q108.EP;q97=q98.covariance;q94=q97.response.inputs.source
    class Phase:
        def normal(self,j):return E(int(j==0))
        def matrix(self,j):return E(2*int(j==0))
        def moving_normal(self,j):return E(3*int(j==0))
    target=post.current.current.state.target
    wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
    cov=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
    pole=post.current.fast.rows.pole+2
    pm=q98.ShiftedPressure(cov,c0=2,c1=4,shape=[1]*13,time_inverse=lambda n,p:1,maximum=11,pole=pole)
    total=q101.FinalRows(q100.SourcedRowChanges(q99.ChangedScalars(pm)),maximum=10,pole=pole,c0=2,c1=4,shape=[1]*13,time_inverse=lambda n,p:1)
    fresh=q102.PostfastSlow(total);reset=q103.CurrentPostslow(fresh)
    final=q104.FinalRadial(reset,maximum=1,pole=pole,c0=2,c1=4,shape=[1]*3,time_inverse=lambda n,p:1)
    rows=q105.ProjectedRows(final)
    return rows,pm,pole

class Tests(unittest.TestCase):
    def test_full_angular_identity_keeps_haar_and_full_differences_once(self):
        y=s.symbols('y',real=True);old=3+s.cos(y);delta=5+s.sin(2*y)
        pre=old+delta;haar=8;RF=7+s.cos(3*y);RS=11+s.sin(y)
        final=s.expand(pre-(pre-haar)+RF+RS)
        self.assertEqual(s.simplify(final-(8+RF+RS)),0)
        self.assertNotEqual(s.simplify(final-(RF+RS)),0)
        self.assertNotEqual(s.simplify(final-(8+2*(RF+RS))),0)
    def test_independent_unweighted_join_and_rate_minimum(self):
        E=q108.EP
        self.assertEqual(q108.join(2,3,5,7,11,13).record(),E(2+33+13*12).record())
        for k in (Q(1,100000),Q(1,5),Q(249,1000)):
            r=q108.rates(k);self.assertEqual(r['full_angular'],Q(3,2)-k)
            self.assertGreater(r['old']-r['sourced'],0)
            self.assertGreater(r['fast']-r['sourced'],0)
            self.assertGreater(r['slow']-r['sourced'],0)
        for k in (0,Q(1,4)):
            with self.assertRaises(ValueError):q108.rates(k)
    def test_same_factory_incoming_program_and_actual_changed_pressure(self):
        rows,pm,pole=setup();a=q108.FullAngular(rows);E=q108.EP;calls=[]
        c=a.means.current;b=c.initial.budget;state=a.state
        vector=lambda v:dict(x=E(v),u=E(v+1),z=E(v+2))
        c.old_full=lambda n:(calls.append(('old',n)) or vector(1))
        c.fast_unweighted=lambda n:(calls.append(('first_fast',n)) or vector(4))
        c.slow=lambda n,w:(calls.append(('first_slow',n,w)) or vector(7))
        b._data=lambda n,p:{k:E(1) for k in ('b','V','G','graph','r')}
        state.target.source.haar.rows=lambda n:(calls.append(('initial_haar',n)) or dict(Etheta=E(2)))
        original_defect=q108.native.covariance.defect;original_edge=q108.radial.edge
        def defect(st,n,**kw):
            self.assertIs(st,state);self.assertEqual(kw,dict(target_rate=Q(2)));calls.append(('signed',n))
            return dict(order=n,pole=2,signed=E(3))
        q108.native.covariance.defect=defect;q108.radial.edge=lambda p,w:E(p+1)
        # Preserve the real Q98 stage/pole/order checks, replace operator leaves only.
        pm.operator=lambda n,p:(calls.append(('new_pressure',n,p)) or dict(KT=E(1),KA=E(1)))
        bundle=dict(source_state=state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=12,covariance_pole=pole,change=[E(1)]*13)
        rows.bound=lambda n,b:dict(source_state=state,order=n,pole=pole+1,fast_difference=dict(Etheta=E(5)),slow_difference=dict(Etheta=E(7)))
        try:
            old=a.old(0)
            self.assertEqual(calls[:5],[('old',3),('first_fast',3),('first_slow',3,False),('initial_haar',0),('signed',1)])
            fixed={k:E(1) for k in ('b','V','G','graph','r')}
            first=q108.angular.angular_difference(0,vector(1),vector(4),fixed)
            second=q108.angular.angular_difference(0,{k:vector(1)[k]+vector(4)[k] for k in vector(1)},vector(7),fixed)
            initial=E(2)*(2*state.target.source.haar.pole+1)
            H=initial+4*E(3)*3
            self.assertEqual(old['coefficient'].record(),(243*(H+first+second)).record())
            out=a.final(0,bundle);pre=out['pre_fast']
            self.assertEqual(out['coefficient'].record(),(pre['old']['coefficient']+(pole+1)*q108.cv(pre['sourced_haar']['theta'])+(pole+2)*E(12)).record())
            self.assertIn(('new_pressure',1,pole),calls)
            self.assertTrue(out['full_angular_row']);self.assertFalse(out['full_axial_row_supplied'])
            for mutation in (dict(order=11),dict(means_updated=True),dict(source_state=object()),dict(stage='postslow')):
                with self.assertRaises(ValueError):a.final(0,dict(bundle,**mutation))
            rows.bound=lambda n,b:dict(source_state=object(),order=n,pole=pole)
            with self.assertRaises(ValueError):a.final(0,bundle)
            a.state=object()
            with self.assertRaises(ValueError):a.old(0)
        finally:q108.native.covariance.defect=original_defect;q108.radial.edge=original_edge
        with self.assertRaises(TypeError):q108.FullAngular(None)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
