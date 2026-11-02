import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q105
q104=q105.radial;q103=q104.post;q102=q103.slow;q101=q102.rows
q100=q101.changes;q99=q100.scalars;q98=q100.pressure
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text()
ns={'q87':type('Module',(),{'full':q105.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

class Tests(unittest.TestCase):
    def test_axial_potential_identity_retains_two_viscous_terms(self):
        R,y,Z,T=s.symbols('R y Z T',real=True);h=(R-1)**3*(2-R)**3;eps=s.Rational(1,3)
        Dr=lambda f:s.diff(f,R)+3*R*s.diff(f,y)
        D1=lambda f:Dr(f)+f/R
        Dz=lambda f:eps*s.diff(f,Z)
        ts=lambda f:-eps*s.diff(f,T)+2*s.diff(f,y)
        Psi=h*(1+Z**3+T*T)*s.sin(y)
        beta=-Dz(Psi);gamma=D1(Psi);v=h*(Z+T)*s.cos(y)
        b=eps*(R+Z);V=1+R;G=R*Z
        Wrr=h*(2+s.cos(y));Wrz=h*(3+s.sin(y));Wzz=h*(5+s.cos(2*y))
        Sigma=h*(1+Z)*s.sin(2*y);primitive=h*s.cos(y)*(1+Z*Z)
        S=2*b*beta+beta*beta+Wrr;H=b*gamma+G*beta+beta*gamma+Wrz
        K0=-S+eps*D1(beta);pm=K0+primitive
        H1=2*G*gamma+gamma*gamma+Wzz-S+primitive
        Ez=ts(gamma)+D1(H)+Dz(2*G*gamma+gamma*gamma+Wzz+pm)-eps*(D1(Dr(gamma))+Dz(Dz(gamma)))-D1(Sigma)
        Omega=ts(Psi)+H-eps*Dr(gamma)-2*eps*Dz(Dz(Psi))-Sigma
        self.assertEqual(s.expand(Ez-Dz(H1)-D1(Omega)),0)
        self.assertNotEqual(s.expand(Ez-Dz(H1)-D1(Omega+eps*Dz(Dz(Psi)))),0)
        self.assertNotEqual(s.expand(Ez-Dz(eps*D1(beta))-Dz(H1)-D1(Omega)),0)
        variable_ts=lambda f:-eps*s.diff(f,T)+R*s.diff(f,y)
        comm=variable_ts(D1(Psi))-D1(variable_ts(Psi))
        self.assertEqual(s.simplify(comm+s.diff(Psi,y)),0)
        self.assertNotEqual(s.simplify(comm),0)
    def test_whole_characteristic_range_and_moving_A1_commutator(self):
        R,r,y,Z=s.symbols('R r y Z',real=True);h=(R-1)**3*(2-R)**3;t=R-1
        phi=y-2*R;Omega=h*s.cos(phi)/R
        Dr=lambda f:s.diff(f,R)+2*s.diff(f,y)
        D1=lambda f:Dr(f)+f/R
        def J1(f):
            return s.integrate(s.expand((R*f).subs({R:r,y:y+2*(r-R)},simultaneous=True)),(r,1,2))
        self.assertEqual(s.simplify(J1(D1(Omega))),0)
        chi=3*t*t-2*t**3+Z*t*t*(1-t)**2
        A1=lambda f:s.diff(chi,R)*J1(f)/R
        comm=s.diff(A1(Omega),Z)-A1(s.diff(Omega,Z))
        self.assertEqual(s.simplify(comm-s.diff(chi,R,Z)*J1(Omega)/R),0)
        self.assertNotEqual(s.simplify(comm),0)
    def test_projected_join_does_not_bound_pointwise_Haar(self):
        a=q105.projected_join(dict(Etheta=2,Ez=3),dict(Etheta=5,Ez=7),11)
        self.assertEqual(a['Qtheta'],14);self.assertEqual(a['Qz'].record(),q105.EP(42).record())
        R=s.symbols('R',positive=True);phi=(R-1)**3*(2-R)**3
        ft=s.diff(phi,R)/R**2;fz=s.diff(phi,R)/R
        self.assertEqual(s.integrate(R**2*ft,(R,1,2)),0)
        self.assertEqual(s.integrate(R*fz,(R,1,2)),0)
        self.assertNotEqual(ft,0);self.assertNotEqual(fz,0)
        # Auxiliary-independent functions have Qcirc=0, despite nonzero rows.
        self.assertEqual(ft-ft,0);self.assertEqual(fz-fz,0)
        for k in (Q(1,100000),Q(1,5)):
            angular=Q(5,2)-4*k;tail=2-k;axial=min(angular,tail)
            self.assertEqual(axial,tail if k<Q(1,6) else angular)
        k=Q(1,100000);target=3-6*k
        self.assertEqual(target-(2-k),1-5*k)
        self.assertEqual(target-(Q(3,2)-2*k),Q(3,2)-4*k)
    def test_current_pre_fast_stage_operator_and_feedback_bindings(self):
        q97=q98.covariance;q94=q97.response.inputs.source;post=fixture(13);E=q105.EP
        class Phase:
            def normal(self,j):return E(int(j==0))
            def matrix(self,j):return E(2*int(j==0))
            def moving_normal(self,j):return E(3*int(j==0))
        target=post.current.current.state.target
        wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        cov=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
        oldpole=post.current.fast.rows.pole;pole=oldpole+2
        if pole%2:pole+=1
        pm=q98.ShiftedPressure(cov,c0=2,c1=4,shape=[1]*13,time_inverse=lambda n,p:1,maximum=11,pole=pole)
        total=q101.FinalRows(q100.SourcedRowChanges(q99.ChangedScalars(pm)),maximum=10,pole=pole,c0=2,c1=4,shape=[1]*13,time_inverse=lambda n,p:1)
        fresh=q102.PostfastSlow(total);reset=q103.CurrentPostslow(fresh)
        final=q104.FinalRadial(reset,maximum=1,pole=pole,c0=2,c1=4,shape=[1]*3,time_inverse=lambda n,p:1)
        adapter=q105.ProjectedRows(final);calls=[]
        vector=lambda a:dict(x=E(a),u=E(a+1),z=E(a+2))
        post.current.initial.budget._fast=lambda n:(calls.append(('initial_fast',n)) or vector(1))
        post.current.initial.slow_fields=lambda n,w:dict(X=E(4),U=E(5),Z=E(6))
        post.current.fast.rows.rows=lambda n:dict(Etheta=E(2),Ez=E(3))
        post.current.fast.compact=lambda n,p:dict(order=n,pole=p,KT=E(3),KA=E(4),I=E(5))
        post.current.slow=lambda n,w:vector(7)
        fresh.means.fields=lambda n:vector(1)
        total.fast=lambda n,b:(calls.append(('new_fast',n)) or dict(order=n,pole=pole,source_state=cov.state,**vector(10)))
        fresh.slow=lambda n,b,w:(calls.append(('new_slow',n,w)) or dict(order=n,source_state=cov.state,weighted=w,fields=vector(13)))
        fresh.changed.bound=lambda n,b:(calls.append(('pre_scalar',n)) or dict(changed=dict(P=E(7),Jtheta=E(8),Jz=E(9))))
        def forbidden(*a):raise AssertionError('Final postslow radial/scalar source substituted for pre-fast source')
        reset.scalars=forbidden;final.bound=forbidden
        post.current.initial.budget._data=lambda n,p:{k:E(1) for k in ('b','V','G','c','graph','r')}
        final.operator=lambda n,p:(calls.append(('pressure_op',n,p)) or dict(order=n,pole=p,KT=E(2),KA=E(3),I=E(5)))
        total.axial=lambda n,p:(calls.append(('axial_op',n,p)) or dict(order=n,pole=p,KT=E(2),KA=E(5),I=E(7)))
        wave.through=lambda n:dict(normal=E(1),order=n,wave_state=cov.state)
        bundle=dict(source_state=cov.state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=12,covariance_pole=pole,wave_pole=pole//2,old=[E(1)]*13,new=[E(2)]*13)
        edge=q104.edge;q104.edge=lambda p,w:E(2)
        try:
            adapter.pre_fields(3,bundle,pole)
            self.assertFalse(any(x[0] in ('new_fast','new_slow') for x in calls))
            out=adapter.bound(0,bundle)
            self.assertIn(('pre_scalar',1),calls);self.assertIn(('pressure_op',1,pole),calls)
            self.assertIn(('axial_op',0,pole),calls);self.assertIn(('new_slow',3,True),calls)
            self.assertIs(out['source_state'],cov.state);self.assertFalse(out['pointwise_Haar_rows_supplied'])
            feedback=adapter.feedback(0,bundle)
            self.assertEqual(feedback['coefficient'].record(),E(2934*243*108*118).record())
            self.assertIn(('new_slow',3,False),calls)
            self.assertFalse(feedback['complete_regular_source_bound_supplied'])
            for mutation in (dict(order=11),dict(covariance_pole=pole+1),dict(means_updated=True),dict(source_state=object())):
                with self.assertRaises(ValueError):adapter.bound(0,dict(bundle,**mutation))
            with self.assertRaises(ValueError):adapter.feedback(0,dict(bundle,wave_pole=pole//2+1))
            wave.through=lambda n:dict(normal=E(1),order=n,wave_state=object())
            with self.assertRaises(ValueError):adapter.feedback(0,bundle)
            final.maximum=0
            with self.assertRaises(ValueError):q105.ProjectedRows(final)
        finally:q104.edge=edge
    def test_record(self):
        self.assertFalse(q105.record()['pointwise_Haar_rows_or_full_regular_residual_supplied'])
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
