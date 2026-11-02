import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q104
q103=q104.post;q102=q103.slow;q101=q102.rows;q100=q101.changes;q99=q100.scalars;q98=q100.pressure
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q104.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

class Tests(unittest.TestCase):
    def test_independent_radial_split_keeps_full_vector_viscosity(self):
        R,y,Z,T=s.symbols('R y Z T',real=True);h=(R-1)**3*(2-R)**3;eps=s.Rational(1,3)
        Dr=lambda f:s.diff(f,R)+2*R*s.diff(f,y)
        Dz=lambda f:eps*s.diff(f,Z)
        D1=lambda f:Dr(f)+f/R
        ts=lambda f:-eps*s.diff(f,T)+2*s.diff(f,y)
        beta=h*(1+Z)*s.cos(y);v=h*(1+T)*s.sin(y);gamma=h*(1+Z*Z)*s.cos(2*y)
        b=R*Z;V=1+R;G=R+T
        Wrr=h*(2+s.cos(y));Wrz=h*(3+s.sin(2*y));Wtt=h*(4+s.cos(2*y))
        S=2*b*beta+beta*beta+Wrr;H=b*gamma+G*beta+beta*gamma+Wrz;C=2*V*v+v*v+Wtt
        vector_vis=Dr(Dr(beta))+Dr(beta)/R-beta/R**2
        gr=-ts(beta)-D1(S)-Dz(H)+C/R+eps*(vector_vis+Dz(Dz(beta)))
        r0=-ts(beta)-S/R-Dz(H)+C/R+eps*Dz(Dz(beta))
        K=-S+eps*D1(beta)
        self.assertEqual(s.simplify(s.expand(gr-r0-Dr(K))),0)
        self.assertNotEqual(s.simplify(gr-r0-Dr(-S)),0)
        self.assertNotEqual(s.simplify(gr-r0-Dr(K+Wrr)),0)
        self.assertEqual(s.simplify(vector_vis-(D1(Dr(beta))-beta/R**2)),0)
    def test_shifted_pressure_and_surviving_A0_at_zero_Haar_moment(self):
        R,r,y,Z=s.symbols('R r y Z',real=True);h=(R-1)**3*(2-R)**3;t=R-1
        chi=3*t*t-2*t**3+Z*t*t*(1-t)**2;rho=30*t*t*(1-t)**2
        K=h*s.cos(y-2*R);r0=h*(s.cos(y-2*R)+Z)
        Dr=lambda f:s.diff(f,R)+2*s.diff(f,y)
        def integral(f,a,b):
            return s.integrate(s.expand(f.subs({R:r,y:y+2*(r-R)},simultaneous=True)),(r,a,b))
        whole=lambda f:integral(f,1,2)
        T0=lambda f:integral(f,1,R)-chi*whole(f)
        A0=lambda f:s.diff(chi,R)*whole(f)
        P=Z*s.integrate(h,(R,1,2));gr=r0+Dr(K)
        pm=K+T0(r0-rho*P)
        self.assertEqual(s.simplify(T0(Dr(K))-K),0)
        self.assertEqual(s.simplify(A0(Dr(K))),0)
        self.assertEqual(s.simplify(pm-T0(gr-rho*P)),0)
        self.assertEqual(s.simplify(Dr(pm)-gr+rho*P+A0(r0-rho*P)),0)
        haar_residual=s.integrate(s.expand(Dr(pm)-gr),(y,0,2*s.pi))/(2*s.pi)
        self.assertEqual(s.simplify(haar_residual+rho*P),0)
        self.assertEqual(s.simplify(s.integrate(haar_residual,(R,1,2))+P),0)
        self.assertNotEqual(s.simplify(pm-T0(r0-rho*P)),0)
        self.assertNotEqual(s.simplify(A0(r0-rho*P).subs(Z,0)),0)
    def test_total_wave_covariance_six_mean_sum_and_positive_envelope(self):
        E=q104.EP;b=dict(old=[E(1),E(2),E(3)],new=[E(4),E(5),E(6)])
        self.assertEqual(q104.total_covariance(2,b).record(),E(2725).record())
        self.assertNotEqual(q104.total_covariance(0,b).record(),E(24).record())
        v=q104.field_sum([('a',0,dict(x=1,u=2,z=3)),('b',2,dict(x=4,u=5,z=6))],2,3)
        for k,a in dict(x=16,u=26,z=36).items():self.assertEqual(v[k].record(),E(243*a).record())
        a=q104.radial_envelope(0,dict(x=1,u=2,z=3),4,1,dict(b=1,V=2,G=3,c=4,graph=5),2)
        for k,n in dict(S=12,H=28,C=36,local_pressure=36,reduced_source=114,direct_source=754,radial_viscosity=512).items():self.assertEqual(a[k].record(),E(n).record())
        with self.assertRaises(ValueError):q104.total_covariance(3,b)
        k=Q(1,100000);target=3-6*k
        self.assertEqual(target-(1-k),2-5*k)
        self.assertEqual(target-(Q(5,2)-k),Q(1,2)-5*k)
        self.assertGreater(target-(Q(5,2)-k),0)
    def test_typed_complete_final_state_and_operator_domains(self):
        q97=q98.covariance;q94=q97.response.inputs.source;post=fixture(12);E=q104.EP
        class Phase:
            def normal(self,j):return E(int(j==0))
            def matrix(self,j):return E(2*int(j==0))
            def moving_normal(self,j):return E(3*int(j==0))
        target=post.current.current.state.target
        wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        cov=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
        oldpole=post.current.fast.rows.pole;pole=oldpole+1
        pm=q98.ShiftedPressure(cov,c0=2,c1=4,shape=[1]*12,time_inverse=lambda n,p:1,maximum=10,pole=pole)
        total=q101.FinalRows(q100.SourcedRowChanges(q99.ChangedScalars(pm)),maximum=9,pole=pole,c0=2,c1=4,shape=[1]*12,time_inverse=lambda n,p:1)
        fresh=q102.PostfastSlow(total);reset=q103.CurrentPostslow(fresh)
        adapter=q104.FinalRadial(reset,maximum=0,pole=pole,c0=2,c1=4,shape=[1]*2,time_inverse=lambda n,p:1)
        calls=[];vector=lambda a:dict(x=E(a),u=E(a+1),z=E(a+2))
        # Independent manufactured lower leaves; no canonical norm certificate.
        post.current.initial.budget._fast=lambda n:(calls.append(('initial_fast',n)) or vector(1))
        post.current.initial.slow_fields=lambda n,w:(calls.append(('initial_slow',n,w)) or dict(X=E(4),U=E(5),Z=E(6)))
        post.current.fast.rows.rows=lambda n:dict(Etheta=E(2),Ez=E(3))
        post.current.fast.compact=lambda n,p:dict(order=n,pole=p,KT=E(3),KA=E(4),I=E(5))
        post.current.slow=lambda n,w:(calls.append(('first_slow',n,w)) or vector(7))
        total.fast=lambda n,b:(calls.append(('new_fast',n)) or dict(order=n,pole=pole,source_state=cov.state,**vector(10)))
        fresh.slow=lambda n,b,w:(calls.append(('new_slow',n,w)) or dict(order=n,source_state=cov.state,weighted=w,fields=vector(13)))
        reset.scalars=lambda n,b:dict(scalars=dict(P=E(1),Jtheta=E(2),Jz=E(3)))
        post.current.initial.budget._data=lambda n,p:{k:E(1) for k in ('b','V','G','c','graph')}
        bundle=dict(source_state=cov.state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=11,covariance_pole=pole,old=[E(1)]*12,new=[E(2)]*12)
        out=adapter.bound(0,bundle)
        self.assertEqual(out['pole'],pole)
        self.assertEqual(out['total_wave_covariance'].record(),E(549).record())
        self.assertEqual({p[0] for p in out['mean']['parts']},{'initial_fast','initial_slow','first_real_fast','first_real_slow','sourced_fast','sourced_slow'})
        self.assertIn(('new_slow',2,True),calls);self.assertIn(('initial_slow',2,True),calls)
        self.assertIn(('initial_fast',2),calls);self.assertFalse(out['complete_residual_improvement'])
        for mutation in (dict(order=10),dict(covariance_pole=pole+1),dict(means_updated=True),dict(source_state=object())):
            with self.assertRaises(ValueError):adapter.bound(0,dict(bundle,**mutation))
        with self.assertRaises(TypeError):q104.FinalRadial(None,maximum=0,pole=1,c0=2,c1=4,shape=[1]*2,time_inverse=lambda n,p:1)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
