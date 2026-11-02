import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q99
q98=q99.pressure
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q99.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

class Tests(unittest.TestCase):
    def test_complete_integrated_rows_with_pressure(self):
        R,z,y=s.symbols('R z y',real=True);h=(R-1)**3*(3-R)**3
        rr=h*(1+z);tt=h*(2+z*z);rz=h*z;zz=h*(3+z);tz=h*(1+z*z)
        gr=-s.diff(rr,R)-rr/R-s.diff(rz,z)+tt/R
        integ=lambda f:s.integrate(f,(R,1,3))
        P=integ(gr);rho=s.Rational(1,2)
        pm=s.integrate(gr-rho*P,R);pm=pm-pm.subs(R,1)
        Et=s.diff(h*(z+1),R)+2*h*(z+1)/R+s.diff(tz,z)
        Ez=s.diff(rz,R)+rz/R+s.diff(zz+pm,z)
        Jt=integ(R*R*tz);Jz=integ(R*zz)-integ(R*R*gr)/2
        c=integ(R*R*rho)/2
        self.assertEqual(s.simplify(integ(R*R*Et)-s.diff(Jt,z)),0)
        self.assertEqual(s.simplify(integ(R*Ez)-s.diff(Jz+c*P,z)),0)
        moving_c=1+z*z;moving_P=z+z*z
        self.assertNotEqual(s.diff(moving_c*moving_P,z),moving_c*s.diff(moving_P,z))
        # The nonzero Haar row is invisible to all three integrated scalar changes.
        mode=(s.diff(h,R)+2*h/R)*s.cos(y)-7*h*s.sin(y)
        self.assertEqual(s.integrate(mode,(y,0,2*s.pi)),0)
        self.assertNotEqual(s.simplify(mode.subs(y,0)),0)
    def test_positive_join_keeps_old_state_and_whole_product(self):
        E=q99.EP
        moments={k:[E(a),E(b)] for k,a,b in [('P',1,2),('Jtheta',3,4),('Jz',5,6)]}
        moments['whole_axial']=[E(62)]
        v=q99.join(0,dict(P=10,Jtheta=20,Jz=30),40,moments)
        for k,a in dict(P=21,Jtheta=43,Jz=65).items():self.assertEqual(v['changed'][k].record(),E(a).record())
        self.assertEqual(v['integrated_theta'].record(),E(24).record())
        self.assertEqual(v['integrated_axial'].record(),E(102).record())
        self.assertNotEqual(v['changed']['P'].record(),v['delta']['P'].record())
        with self.assertRaises(ValueError):q99.join(1,dict(P=1,Jtheta=1,Jz=1),1,moments)
    def test_typed_current_postslow_not_stale_signed_scalars(self):
        q97=q98.covariance;q94=q97.response.inputs.source;post=fixture(3)
        class Phase:
            def normal(self,j):return q99.EP(int(j==0))
            def matrix(self,j):return q99.EP(2*int(j==0))
            def moving_normal(self,j):return q99.EP(3*int(j==0))
        target=post.current.current.state.target
        wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        changed=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
        pm=q98.ShiftedPressure(changed,c0=2,c1=4,shape=[1]*4,time_inverse=lambda n,p:1,maximum=2,pole=40)
        adapter=q99.ChangedScalars(pm);calls=[]
        # Manufactured old norm leaves isolate the stage binding, not actual fields.
        post.scalars=lambda n:(calls.append(('postslow',n)) or {k:q99.EP(11) for k in q99.KEYS})
        post.axial=lambda n:(calls.append(('axial',n)) or q99.EP(13))
        bundle=dict(source_state=changed.state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=2,covariance_pole=40,change=[q99.EP(1)]*3)
        v=adapter.bound(0,bundle)
        self.assertEqual(calls,[('postslow',1),('axial',0)])
        self.assertFalse(v['means_updated']);self.assertFalse(v['fresh_mean_source_rows_supplied'])
        self.assertIs(v['source_state'],changed.state)
        with self.assertRaises(ValueError):adapter.bound(0,dict(bundle,order=1))
        with self.assertRaises(ValueError):adapter.bound(0,dict(bundle,means_updated=True))
        with self.assertRaises(TypeError):q99.ChangedScalars(None)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
