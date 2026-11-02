import unittest,json
from types import SimpleNamespace as NS
from fractions import Fraction as Q
import sympy as s
import audit,q107

class Tests(unittest.TestCase):
    def test_fixed_mean_full_row_subtraction_and_no_graph_multiplier(self):
        R,Z,T=s.symbols('R Z T',real=True);e=s.Rational(1,3)
        p=audit.HERE.parent/'core_current_head_q106_pointwise_haar_ledger/q106.py'
        import importlib.util
        spec=importlib.util.spec_from_file_location('haar_ledger',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        means=dict(x=R*Z,u=R*T,z=R*Z**2);prod={k:R**2+Z for k in ('xu','zu','xz','zz')}
        base=dict(b=R,V=R**2,G=Z);stress=dict(t=R*T,z=R*Z)
        W={k:R*(i+Z) for i,k in enumerate(('rr','rt','rz','tt','zt','zz'))}
        D={k:R**2*(i+Z**2) for i,k in enumerate(W)};oldpm=R*Z;dp=R**3*Z
        old=m.rows(R,Z,T,e,means,prod,W,base,oldpm,stress)
        new=m.rows(R,Z,T,e,means,prod,{k:W[k]+D[k] for k in W},base,oldpm+dp,stress)
        theta=s.diff(D['rt'],R)+2*D['rt']/R+e*s.diff(D['zt'],Z)
        axial=s.diff(D['rz'],R)+D['rz']/R+e*s.diff(D['zz']+dp,Z)
        self.assertEqual(s.simplify(new['theta']-old['theta']-theta),0)
        self.assertEqual(s.simplify(new['axial']-old['axial']-axial),0)
        y=s.symbols('y',real=True);g=R**2*(Z+s.cos(y));Dr=lambda f:s.diff(f,R)+17*R*s.diff(f,y)
        avg=lambda f:s.integrate(s.expand(f),(y,-s.pi,s.pi))/(2*s.pi)
        self.assertEqual(s.simplify(avg(Dr(g))-s.diff(avg(g),R)),0)
    def test_compensated_haar_pressure_zero_integral_and_moving_rho(self):
        R,r,Z=s.symbols('R r Z',real=True);e=s.Rational(1,3);t=R-1;h=t**3*(1-t)**3
        rho=30*t*t*(1-t)**2+Z*s.diff(h,R)
        whole=lambda f:s.integrate(s.expand(f).subs(R,r),(r,1,2))
        left=lambda f:s.integrate(s.expand(f).subs(R,r),(r,1,R))
        rr=R*h*(1+Z*Z);rz=h*Z;tt=2*R*h*(1+Z)
        H=(tt-rr)/R-e*s.diff(rz,Z);P=whole(H);f=H-rho*P
        self.assertEqual(s.simplify(whole(rho)),1);self.assertEqual(s.simplify(whole(f)),0)
        chi=3*t*t-2*t**3+Z*t*t*(1-t)**2
        pm=-rr+left(f)-chi*whole(f)
        gr=-s.diff(rr,R)-rr/R-e*s.diff(rz,Z)+tt/R
        self.assertEqual(s.simplify(s.diff(pm,R)-gr+rho*P),0)
        self.assertEqual(s.simplify(pm+rr-left(f)),0)
        self.assertNotEqual(s.simplify(s.diff(pm,Z)-(-s.diff(rr,Z)+left(s.diff(H,Z)-rho*s.diff(P,Z)))),0)
        self.assertNotEqual(s.simplify(s.diff(left(f),R)-gr+rho*P),0)
    def test_opposite_modes_cross_and_square_survive_haar(self):
        y,l=s.symbols('y l',real=True)
        old=s.Matrix([s.cos(y),s.cos(y),s.sin(2*y)])
        new=s.Matrix([l*s.cos(y),2*l*s.cos(y),l*s.sin(2*y)])
        avg=lambda f:s.integrate(s.expand_trig(f),(y,-s.pi,s.pi))/(2*s.pi)
        delta=((old+new)*(old+new).T-old*old.T).applyfunc(avg)
        cross=(old*new.T+new*old.T).applyfunc(avg);square=(new*new.T).applyfunc(avg)
        self.assertEqual(s.simplify(delta-cross-square),s.zeros(3))
        self.assertEqual(s.expand(delta[0,1]),3*l/2+l*l)
        self.assertEqual(avg(new[0]),0);self.assertNotEqual(delta[0,1],0)
        self.assertNotEqual(s.simplify(delta-square),s.zeros(3))
    def test_exact_envelope_and_grades(self):
        E=q107.EP;a=q107.envelope(0,[1,2],1,3)
        # Through-0=1, through-1=1+5*2=11.
        for k,v in dict(theta=24,axial=26,radial_flux=13,slow_covariance=11,slow_pressure=3).items():self.assertEqual(a[k].record(),E(v).record())
        for k in (Q(1,100000),Q(1,5)):
            r=q107.rates(k);self.assertEqual(r['changed_haar']-r['full_changed_row'],k)
            self.assertGreater(r['remaining_target_gap'],0)
        for k in (0,Q(1,4)):
            with self.assertRaises(ValueError):q107.rates(k)
        with self.assertRaises(ValueError):q107.envelope(1,[1,2],1,3)
    def test_actual_q98_protocol_and_same_state_guards(self):
        # Manufactured operator leaves; real Q98 bound/envelope/guards execute.
        q98=q107.pressure;state=NS(target=NS(width=Q(1,4)));cp=NS(Rmin=Q(1),Rmax=Q(2),weighted_rho=lambda n:1)
        pm=object.__new__(q98.ShiftedPressure);pm.state=state;pm.maps=NS(profile=cp)
        pm.operator=lambda n,p:dict(KT=1,KA=1)
        c=object.__new__(q107.changes.SourcedRowChanges);c.state=state;c.pressure=pm
        a=q107.HaarChanges(c)
        bundle=dict(source_state=state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=2,covariance_pole=2,change=[1,2,3])
        out=a.bound(0,bundle);self.assertIs(out['source_state'],state)
        self.assertEqual(out['pressure_request'],1);self.assertEqual(out['covariance_request'],2)
        self.assertFalse(out['graph_multiplier_required']);self.assertTrue(out['incoming_haar_rows_still_required'])
        for mutation in (dict(source_state=object()),dict(means_updated=True),dict(order=1),dict(stage='postslow')):
            with self.assertRaises(ValueError):a.bound(0,dict(bundle,**mutation))
        pm.state=object()
        with self.assertRaises(ValueError):a.bound(0,bundle)
        with self.assertRaises(TypeError):q107.HaarChanges(None)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
