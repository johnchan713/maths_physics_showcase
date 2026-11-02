import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q103
q102=q103.slow;q101=q102.rows;q100=q101.changes;q99=q100.scalars;q98=q100.pressure
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text();ns={'q87':type('Module',(),{'full':q103.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

class Tests(unittest.TestCase):
    def test_current_haar_keeps_both_slow_histories(self):
        y,Z=s.symbols('y Z',real=True);initial=1+Z;previous=2+Z*Z
        old=initial+previous+3*s.cos(y)+5*s.sin(2*y);inc=3+Z
        haar=lambda f:s.integrate(f,(y,0,2*s.pi))/(2*s.pi)
        self.assertEqual(s.simplify(haar(old)-initial-previous),0)
        self.assertEqual(s.simplify(haar(old*inc)-(initial+previous)*inc),0)
        self.assertNotEqual(s.simplify(haar(old*inc)-initial*inc),0)
        v=q103.haar_sum(dict(X=1,U=2,Z=3),dict(x=4,u=5,z=6))
        for k,a in dict(x=5,u=7,z=9).items():self.assertEqual(v[k].record(),q103.EP(243*a).record())
    def test_complete_remainder_polynomial(self):
        E=q103.EP;old=dict(x=E(1),u=E(2),z=E(3));inc=dict(x=E(5),u=E(7),z=E(11))
        v=q103.reset(0,old,inc,E(2),3)
        for k,a in dict(P=750,Jtheta=360,Jz=936).items():self.assertEqual(v[k].record(),E(a).record())
        stale=q103.reset(0,{k:E() for k in old},inc,E(2),3)
        self.assertNotEqual(v['Jtheta'].record(),stale['Jtheta'].record())
    def test_whole_potential_G0_is_an_actual_gate(self):
        R=s.symbols('R',real=True);phi=(R-1)**3*(3-R)**3
        # Preserve the zero angular velocity moment as well as compact support.
        c=s.integrate(R**3*phi,(R,1,3))/s.integrate(R**2*phi,(R,1,3))
        u=phi*(R-c)
        self.assertEqual(s.integrate(R*R*u,(R,1,3)),0)
        # G(R)=R leaves a nonzero leak; the G=0 remainder K is still zero.
        leak=s.integrate(R**3*u,(R,1,3))
        self.assertGreater(leak,0)
        v=q103.reset(0,{k:q103.EP() for k in ('x','u','z')},dict(x=q103.EP(),u=q103.EP(1),z=q103.EP()),q103.EP(1),1)
        self.assertEqual(v['Jtheta'].record(),q103.EP().record())
    def test_radial_viscosity_cancels_only_in_compact_scalar_integrals(self):
        R=s.symbols('R',positive=True);x=(R-1)**3*(3-R)**3
        visc=s.diff(s.diff(x,R)+x/R,R)
        self.assertNotEqual(s.simplify(visc),0)
        self.assertEqual(s.integrate(visc,(R,1,3)),0)
        self.assertEqual(s.integrate(R*R*visc,(R,1,3)),0)
        noncompact=s.diff(s.diff(R*R,R)+R,R)
        self.assertNotEqual(s.integrate(noncompact,(R,1,3)),0)
    def test_typed_current_history_new_slow_and_pressure_products(self):
        q97=q98.covariance;q94=q97.response.inputs.source;post=fixture(13);E=q103.EP
        class Phase:
            def normal(self,j):return E(int(j==0))
            def matrix(self,j):return E(2*int(j==0))
            def moving_normal(self,j):return E(3*int(j==0))
        target=post.current.current.state.target
        wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        cov=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
        pole=post.current.fast.rows.pole+1
        pm=q98.ShiftedPressure(cov,c0=2,c1=4,shape=[1]*13,time_inverse=lambda n,p:1,maximum=11,pole=pole)
        total=q101.FinalRows(q100.SourcedRowChanges(q99.ChangedScalars(pm)),maximum=10,pole=pole,c0=2,c1=4,shape=[1]*13,time_inverse=lambda n,p:1)
        fresh=q102.PostfastSlow(total);adapter=q103.CurrentPostslow(fresh);calls=[]
        # Distinct manufactured leaves expose an accidental first-real reuse.
        post.current.initial.slow_fields=lambda n,w:(calls.append(('initial',n,w)) or dict(X=E(1),U=E(2),Z=E(3)))
        post.current.slow=lambda n,w:(calls.append(('previous',n,w)) or dict(x=E(4),u=E(5),z=E(6)))
        fresh.slow=lambda n,b,w:(calls.append(('fresh',n,w)) or dict(fields=dict(x=E(7),u=E(11),z=E(13)),source_state=cov.state,order=n,weighted=w))
        budget=post.current.initial.budget
        budget.fixed=lambda n,p:(calls.append(('base',n,p)) or dict(b=E(2)))
        cp=post.current.current.profile
        cp.c_rho=lambda n:(calls.append(('c_rho',n)) or E(4+n))
        cp.weighted_rho=lambda n:(calls.append(('rho',n)) or E(6+n))
        bundle=dict(source_state=cov.state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=12,covariance_pole=pole,change=[E(1)]*13)
        a=adapter.scalars(0,bundle);an=adapter.scalars(1,bundle)
        self.assertEqual(a['old_haar']['z'].record(),E(243*9).record())
        self.assertIn(('initial',2,False),calls);self.assertIn(('previous',2,False),calls);self.assertIn(('fresh',2,False),calls)
        self.assertIn(('base',2,2*budget.pole),calls)
        self.assertEqual(adapter.axial(0,bundle).record(),(an['scalars']['Jz']+243*5*an['scalars']['P']).record())
        self.assertIn(('c_rho',1),calls)
        self.assertEqual(adapter.compensator(0,bundle).record(),(243*6*a['scalars']['P']).record())
        self.assertFalse(a['actual_native_contracts_certified']);self.assertFalse(a['means_updated'])
        for mutation in (dict(order=10),dict(means_updated=True),dict(source_state=object())):
            with self.assertRaises(ValueError):adapter.scalars(0,dict(bundle,**mutation))
        with self.assertRaises(ValueError):adapter.axial(0,dict(bundle,order=11))
        with self.assertRaises(TypeError):q103.CurrentPostslow(None)
    def test_grade_comparison(self):
        for k in (Q(1,100000),Q(1,5)):
            C=Q(3,2)-k;alpha=2-3*k
            self.assertGreater(alpha,1)
            self.assertEqual(min(C+1,2*C,C+2),C+1)
            self.assertGreater(2*C,C+1)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
