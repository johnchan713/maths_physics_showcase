import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q100
q99=q100.scalars;q98=q100.pressure
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q100.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

class Tests(unittest.TestCase):
    def test_independent_full_rows_keep_pressure_and_connections(self):
        R,y,Z,T=s.symbols('R y Z T',real=True);eps=s.Rational(1,3)
        Dr=lambda f:s.diff(f,R)+3*R*s.diff(f,y)
        Dz=lambda f:eps*s.diff(f,Z)
        De=lambda f,e:Dr(f)+e*f/R
        ts=lambda f:-eps*s.diff(f,T)+2*s.diff(f,y)
        x,u,z=R+T,R*Z+s.cos(y),R*T+s.sin(y)
        b,V,G=R*T,R+Z,2+Z*T
        base=[R+Z,R*T,R*Z];delta=[R*Z*s.cos(y),R*T*s.sin(y),(R+Z)*s.cos(y)]
        pm0=R*R+T;dpm=R*Z*Z*s.sin(y)
        def rows(W,pm):
            rt=(b+x)*(V+u)+W[0];at=(G+z)*(V+u)+W[1]
            az=(b+x)*(G+z)+W[1]
            et=ts(u)+De(rt,2)+Dz(at)-eps*(Dr(Dr(u))+Dr(u)/R-u/R**2+Dz(Dz(u)))
            ez=ts(z)+De(az,1)+Dz((G+z)**2+W[2]+pm)-eps*(De(Dr(z),1)+Dz(Dz(z)))
            return et,ez
        old=rows(base,pm0);new=rows([a+b for a,b in zip(base,delta)],pm0+dpm)
        expected=(De(delta[0],2)+Dz(delta[1]),De(delta[1],1)+Dz(delta[2]+dpm))
        for a,b,c in zip(new,old,expected):self.assertEqual(s.simplify(a-b-c),0)
        self.assertNotEqual(s.simplify(expected[1]-(De(delta[1],1)+Dz(delta[2]))),0)
        # All three new components have zero Haar, yet neither row is zero.
        for f in expected:
            self.assertEqual(s.integrate(f,(y,0,2*s.pi)),0)
            self.assertNotEqual(s.simplify(f),0)
    def test_positive_norm_envelope_and_rate_margin(self):
        v=q100.envelope(0,[1,2],1,3,5)
        for k,a in dict(radial_transport=77,angular_connection=2,axial_connection=1,covariance_axial_gradient=11,pressure_axial_gradient=5,Etheta=90,Ez=94,Qtheta=180,Qz=188).items():
            self.assertEqual(v[k].record(),q100.EP(a).record())
        v=q100.envelope(1,[1,2,3],1,3,5)
        self.assertEqual(v['Etheta'].record(),q100.EP(820).record())
        self.assertEqual(v['Ez'].record(),q100.EP(808).record())
        for k in (Q(1,100000),Q(1,5)):
            self.assertGreater(2-4*k,Q(3,2)-2*k)
        with self.assertRaises(ValueError):q100.envelope(1,[1,2],1,3,5)
    def test_typed_pressure_next_jet_and_stale_fast_refusal(self):
        q97=q98.covariance;q94=q97.response.inputs.source;post=fixture(3)
        class Phase:
            def normal(self,j):return q100.EP(int(j==0))
            def matrix(self,j):return q100.EP(2*int(j==0))
            def moving_normal(self,j):return q100.EP(3*int(j==0))
        target=post.current.current.state.target
        wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        changed=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
        pm=q98.ShiftedPressure(changed,c0=2,c1=4,shape=[1]*4,time_inverse=lambda n,p:1,maximum=2,pole=40)
        adapter=q100.SourcedRowChanges(q99.ChangedScalars(pm));calls=[]
        post.current.initial.budget._data=lambda n,p:(calls.append(('graph',n,p)) or dict(graph=q100.EP(3)))
        original=pm.bound
        def traced(n,b):calls.append(('pressure',n,b['covariance_pole']));return original(n,b)
        pm.bound=traced
        bundle=dict(source_state=changed.state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=2,covariance_pole=40,change=[q100.EP(1)]*3)
        out=adapter.bound(0,bundle)
        self.assertEqual(calls[0],('pressure',1,40));self.assertEqual(calls[1][0:2],('graph',1))
        self.assertEqual(out['Qz'].record(),(2*out['Ez']).record())
        self.assertFalse(out['total_changed_row_bound_supplied'])
        for mutation in (dict(order=1),dict(covariance_pole=41),dict(means_updated=True),dict(source_state=object())):
            with self.assertRaises(ValueError):adapter.bound(0,dict(bundle,**mutation))
        with self.assertRaises(NotImplementedError):adapter.fresh_fast(0,bundle)
        with self.assertRaises(TypeError):q100.SourcedRowChanges(None)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
