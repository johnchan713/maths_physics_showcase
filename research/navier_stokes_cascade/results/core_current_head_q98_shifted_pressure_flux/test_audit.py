import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q98
p=Path(__file__).resolve().parent.parent
text=(p/'core_current_head_q87_final_angular_adapter/test_audit.py').read_text()
ns={'q87':type('Module',(),{'full':q98.full}),'Q':Q}
exec(text[text.index('def fixture('):text.index('class Tests')],ns);fixture=ns['fixture']

class Tests(unittest.TestCase):
    def test_shifted_pressure_split_and_moving_cutoff(self):
        R,t,y,z=s.symbols('R t y z',real=True);I=s.I
        h=(R-1)**3*(2-R)**3
        F=h*s.exp(I*y);H=R+z+(R+z)*s.exp(I*y)
        Dr=lambda f:s.diff(f,R)+2*s.diff(f,y)
        chi=(R-1)+z*(R-1)*(2-R)
        def moment(f,l,u):
            shifted=f.subs({R:t,y:y+2*(t-R)},simultaneous=True)
            return s.integrate(s.expand(shifted),(t,l,u))
        def T(f):return moment(f,1,R)-chi*moment(f,1,2)
        def A(f):return s.diff(chi,R)*moment(f,1,2)
        # These are full characteristic integrals, with nonzero auxiliary mode.
        self.assertEqual(s.simplify(T(Dr(F))-F),0)
        self.assertEqual(s.simplify(A(Dr(F))),0)
        P=s.integrate(R+z,(R,1,2));source=-Dr(F)+H
        direct=T(source-P);split=-F+T(H-P)
        self.assertEqual(s.simplify(direct-split),0)
        self.assertEqual(s.simplify(s.diff(direct-split,z)),0)
        self.assertNotEqual(s.simplify(direct-T(H-P)),0)
        self.assertNotEqual(s.simplify(T(P)),0)
        comm=s.simplify(s.diff(T(H-P),z)-T(s.diff(H-P,z)))
        self.assertNotEqual(comm,0)
        self.assertEqual(s.simplify(comm+s.diff(chi,z)*moment(H-P,1,2)),0)
        residual=s.simplify(Dr(split)-source)
        self.assertEqual(s.simplify(residual+P+A(H-P)),0)
    def test_connection_and_trace_defects_are_not_Haar_shortcuts(self):
        R=s.symbols('R',positive=True);h=(R-1)**3*(2-R)**3
        self.assertEqual(s.integrate(s.diff(h,R),(R,1,2)),0)
        self.assertNotEqual(s.integrate(s.diff(h,R)+h/R,(R,1,2)),0)
        # Opposite nonzero Fourier jumps have zero total jump, yet a cumulative defect.
        y=s.symbols('y',real=True);chi=s.Rational(1,3);jumps=[2*s.cos(y),-2*s.cos(y)]
        A=-sum(jumps);T_defect=-(jumps[0]-chi*sum(jumps))
        self.assertEqual(A,0);self.assertEqual(T_defect,-2*s.cos(y))
        self.assertEqual(s.integrate(T_defect,(y,0,2*s.pi)),0)
        # A piecewise flux with zero exterior traces and jump -1 has
        # classical derivative 1 on a unit hull: delta P gains Haar(J)=-1.
        self.assertEqual(-s.integrate(1,(R,1,2)),-1)
        self.assertEqual(s.integrate(-1+s.cos(y),(y,0,2*s.pi))/(2*s.pi),-1)
    def test_positive_norm_polynomial_and_strict_gain(self):
        v=q98.envelope(0,[1,2],1,2,1,3,5,7)
        for k,n in dict(local_flux=1,H=13,P=13,compensated=52,pressure=261,A0=364,radial_residual=403).items():
            self.assertEqual(v[k].record(),q98.EP(n).record())
        k=Q(1,100000);self.assertGreater(Q(3,2)-k,Q(3,2)-2*k)
        with self.assertRaises(ValueError):q98.envelope(1,[1,2],1,2,1,3,5,7)
    def test_current_binding_exact_pole_and_reserve(self):
        q97=q98.covariance;q94=q97.response.inputs.source;post=fixture(3)
        class Phase:
            def normal(self,j):return q98.EP(int(j==0))
            def matrix(self,j):return q98.EP(2*int(j==0))
            def moving_normal(self,j):return q98.EP(3*int(j==0))
        target=post.current.current.state.target
        wave=q94.waves.CurrentWaves(post,regular_wave=target.root.wave,phase=Phase(),inverse_normal_square=1)
        changed=q97.SourcedCovariance(q97.response.SourcedResponse(q97.response.inputs.PathInputs(q94.CurrentSource(wave))))
        factory=q98.ShiftedPressure(changed,c0=2,c1=4,shape=[1]*4,time_inverse=lambda n,p:1,maximum=2,pole=40)
        # Manufactured norm leaves isolate the operator/state interface, not actual fields.
        bundle=dict(source_state=changed.state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=1,covariance_pole=40,change=[q98.EP(1),q98.EP(2)])
        out=factory.bound(0,bundle)
        self.assertEqual(out['covariance_pole'],40);self.assertEqual(factory.operator.d,post.current.fast.rows.pressure.d)
        self.assertIs(factory.maps.profile,changed.state.target.source.profile)
        self.assertFalse(out['actual_native_contracts_certified'])
        for mutation in (dict(order=0),dict(source_state=object()),dict(means_updated=True),dict(covariance_pole=41)):
            with self.assertRaises(ValueError):factory.bound(0,dict(bundle,**mutation))
        with self.assertRaises(TypeError):q98.ShiftedPressure(None,c0=2,c1=4,shape=[1]*4,time_inverse=lambda n,p:1,maximum=2,pole=40)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
