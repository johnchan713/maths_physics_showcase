import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q102
q101=q102.rows;q100=q101.changes;q99=q100.scalars;q98=q100.pressure
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q102.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

class Tests(unittest.TestCase):
    def test_full_radial_source_and_scalar_correlation_identities(self):
        R,y,Z,T=s.symbols('R y Z T',real=True);h=(R-1)**3*(3-R)**3;eps=s.Rational(1,3)
        Dr=lambda f:s.diff(f,R)+2*s.diff(f,y)
        Dz=lambda f:eps*s.diff(f,Z)
        ts=lambda f:-eps*s.diff(f,T)+3*s.diff(f,y)
        D1=lambda f:Dr(f)+f/R
        old=[h*(a+Z)*s.cos(y) for a in (1,2,3)]
        inc=[h*(a+Z)*s.cos(y) for a in (4,5,6)]
        x,u,z=old;X,U,V=inc;b=R*Z;G=R+Z;B=R+T
        def gr(a):
            x,u,z=a
            return -ts(x)-D1(2*b*x+x*x)-Dz(b*z+G*x+x*z)+(2*B*u+u*u)/R+eps*(Dr(D1(x))+Dz(Dz(x)))
        haar=lambda f:s.integrate(s.expand(f),(y,0,2*s.pi))/(2*s.pi)
        S=haar(2*x*X+X*X);H=haar(x*V+z*X+X*V);A=haar(2*u*U+U*U)
        delta=haar(gr([a+b for a,b in zip(old,inc)])-gr(old))
        self.assertEqual(s.simplify(delta-(-s.diff(S,R)-S/R-Dz(H)+A/R)),0)
        self.assertEqual(s.simplify(delta-((A-S)/R-Dz(H))+s.diff(S,R)),0)
        self.assertEqual(S.subs(R,1),0);self.assertEqual(S.subs(R,3),0)
        K=haar(z*U+u*V+U*V);Zf=haar(2*z*V+V*V)
        flux=haar((G+z+V)*(B+u+U)-(G+z)*(B+u))
        self.assertEqual(s.simplify(flux-K),0)
        fullJ=R*Zf-R*R*delta/2
        reduced=R*Zf-R*(S+A)/2+R*R*Dz(H)/2
        self.assertEqual(s.simplify(fullJ-reduced-s.diff(R*R*S,R)/2),0)
    def test_positive_join_and_moving_whole_pressure_product(self):
        E=q102.EP;old=dict(x=E(2),u=E(3),z=E(5));inc=dict(x=E(7),u=E(11),z=E(13))
        v=q102.correlation_join(0,dict(P=E(10),Jtheta=E(20),Jz=E(30)),old,inc,2)
        for k,a in dict(P=832,Jtheta=474,Jz=1014).items():self.assertEqual(v['fast_change'][k].record(),E(a).record())
        for k,a in dict(P=842,Jtheta=494,Jz=1044).items():self.assertEqual(v['scalars'][k].record(),E(a).record())
        self.assertEqual(q102.whole_axial(0,v['scalars'],4).record(),E(4412).record())
        zero=q102.correlation_join(0,dict(P=E(10),Jtheta=E(20),Jz=E(30)),{k:E() for k in old},inc,2)
        self.assertNotEqual(v['fast_change']['P'].record(),zero['fast_change']['P'].record())
        Z=s.symbols('Z');c=1+Z*Z;P=Z+Z*Z
        self.assertEqual(s.diff(c*P,Z).subs(Z,1),10)
        self.assertEqual((c*s.diff(P,Z)).subs(Z,1),6)
    def test_typed_same_state_new_fast_and_next_scalar_slow_map(self):
        q97=q98.covariance;q94=q97.response.inputs.source;post=fixture(10);E=q102.EP
        class Phase:
            def normal(self,j):return E(int(j==0))
            def matrix(self,j):return E(2*int(j==0))
            def moving_normal(self,j):return E(3*int(j==0))
        target=post.current.current.state.target
        wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        cov=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
        pole=post.current.fast.rows.pole+1
        pm=q98.ShiftedPressure(cov,c0=2,c1=4,shape=[1]*10,time_inverse=lambda n,p:1,maximum=8,pole=pole)
        total=q101.FinalRows(q100.SourcedRowChanges(q99.ChangedScalars(pm)),maximum=7,pole=pole,c0=2,c1=4,shape=[1]*10,time_inverse=lambda n,p:1)
        adapter=q102.PostfastSlow(total);calls=[]
        # Manufactured norm leaves isolate composition; actual native norms remain open.
        total.means.fields=lambda n:(calls.append(('old',n)) or dict(x=E(2),u=E(3),z=E(5)))
        def fast(n,b):
            calls.append(('fast',n));return dict(order=n,pole=pole,source_state=cov.state,x=E(7),u=E(11),z=E(13))
        total.fast=fast
        adapter.changed.bound=lambda n,b:(calls.append(('changed',n)) or dict(changed=dict(P=E(10+n),Jtheta=E(20+n),Jz=E(30+n))))
        def maps(n):return dict(H=E(2 if n==0 else 13),Kp=E(3),Kz=E(5),Kg=E(7),Kpsi=E(11 if n==0 else 17))
        total.current.initial.budget.slow=maps
        bundle=dict(source_state=cov.state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=9,covariance_pole=pole,change=[E(1)]*10)
        previous_edge=q102.edge
        try:
            q102.edge=lambda p,w:(calls.append(('edge',p)) or Q(p+1))
            a=adapter.postfast(0,bundle)['scalars'];an=adapter.postfast(1,bundle)['scalars']
            calls.clear();out=adapter.slow(0,bundle)
            self.assertEqual(out['fields']['u'].record(),(2*(3*a['P']+5*a['Jz'])).record())
            self.assertEqual(out['fields']['z'].record(),(14*a['Jtheta']).record())
            self.assertEqual(out['fields']['x'].record(),(13*17*an['Jtheta']).record())
            self.assertIn(('old',2),calls);self.assertIn(('fast',2),calls);self.assertIn(('changed',1),calls)
            self.assertIn(('edge',pole),calls);self.assertFalse(out['actual_native_contracts_certified'])
            weighted=adapter.slow(0,bundle,True)
            for k,v in out['fields'].items():self.assertEqual(weighted['fields'][k].record(),(q102.cv(total.current.initial.budget.floor)*v).record())
            for mutation in (dict(order=7),dict(source_state=object()),dict(means_updated=True)):
                with self.assertRaises(ValueError):adapter.postfast(0,dict(bundle,**mutation))
            with self.assertRaises(ValueError):adapter.slow(0,dict(bundle,order=8))
            with self.assertRaises(TypeError):adapter.slow(0,bundle,1)
        finally:q102.edge=previous_edge
        with self.assertRaises(TypeError):q102.PostfastSlow(None)
    def test_rate_margins(self):
        for k in (Q(1,100000),Q(1,5)):
            B=1-k;D=Q(3,2)-2*k;C=Q(3,2)-k
            self.assertGreater(B+D,C);self.assertGreater(2*D,C)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
