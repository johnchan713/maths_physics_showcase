import unittest,json
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import audit,q109
text=(audit.HERE.parent/'core_current_head_q108_full_angular_join/test_audit.py').read_text()
ns=dict(q108=q109.angular,Path=Path,Q=Q,__file__=__file__)
exec(text[text.index('q105=q108.projected'):text.index('\nclass Tests(')],ns);setup=ns['setup']

class Tests(unittest.TestCase):
    def test_haar_compact_remainder_zero_but_actual_velocity_correction_nonzero(self):
        R,r,y=s.symbols('R r y',real=True);h=(R-1)**3*(2-R)**3;t=R-1
        chi=3*t*t-2*t**3;f=h*(3+s.cos(y-2*R))
        avg=lambda g:s.integrate(s.expand(g),(y,-s.pi,s.pi))/(2*s.pi)
        whole=lambda g:s.integrate(s.expand((R*g).subs({R:r,y:y+2*(r-R)},simultaneous=True)),(r,1,2))
        A=lambda g:s.diff(chi,R)*whole(g)/R
        Qf=f-avg(f)
        self.assertEqual(s.simplify(avg(A(Qf))),0)
        self.assertNotEqual(s.simplify(A(Qf)),0)
        self.assertNotEqual(s.simplify(A(f)),0)
        self.assertEqual(s.simplify(avg(A(f))-A(avg(f))),0)
        # Removing A1 from the velocity would change its square/correlations.
        self.assertNotEqual(s.simplify(avg((Qf-A(Qf))**2)-avg(Qf**2)),0)
    def test_signed_haar_axial_pressure_local_memory_is_essential(self):
        R,Z=s.symbols('R Z',real=True);h=(R-1)**3*(2-R)**3;e=s.Rational(1,3)
        rr=R*h*Z;rz=h*(1+Z);zz=2*h*Z**2;primitive=h*Z**2
        pm=-rr+primitive
        row=s.diff(rz,R)+rz/R+e*s.diff(zz+pm,Z)
        missing=s.diff(rz,R)+rz/R+e*s.diff(zz+primitive,Z)
        self.assertEqual(s.simplify(row-missing+e*s.diff(rr,Z)),0)
        self.assertNotEqual(s.simplify(row-missing),0)
    def test_grade_minimum_and_raw_jet_padding(self):
        for k in (Q(1,100000),Q(1,5),Q(249,1000)):
            a=q109.rates(k);self.assertEqual(a['full_axial'],Q(3,2)-k)
            self.assertEqual(a['current_A1']-a['full_axial'],Q(1,2))
            self.assertGreater(a['target_gap'],0)
        for k in (0,Q(1,4)):
            with self.assertRaises(ValueError):q109.rates(k)
        # A through-j norm contains 5^j/j! times the raw derivative.
        from math import factorial
        raw=[Q(2),Q(3),Q(5),Q(7)]
        for j in range(4):
            through=sum(Q(5**i,factorial(i))*raw[i] for i in range(j+1))
            self.assertGreaterEqual(factorial(j)*through,raw[j])
    def test_same_state_pressure_pole_requests_and_old_vs_current_A1(self):
        rows,pm,pole=setup();context=q109.angular.FullAngular(rows);a=q109.FullAxial(context)
        E=q109.EP;state=a.state;calls=[];fresh=rows.fresh.rows
        oldbound=fresh.bound;olddefect=q109.angular.native.covariance.defect;oldedge=q109.radial.edge
        fresh.bound=lambda n,b:dict(source_state=state,order=n,pole=pole,fast_difference=dict(Ez=E(5)),slow_difference=dict(Ez=E(7)),post_signed_A1=E(99999))
        fresh.pressure=lambda n,p:(calls.append(('old_pressure',n,p)) or dict(order=n,pole=p,KT=E(2),KA=E(3),I=E(1)))
        def defect(st,n,**kw):
            self.assertIs(st,state);self.assertEqual(kw,dict(target_rate=Q(2)));calls.append(('old_signed',n))
            return dict(order=n,pole=pole-1,signed=E(n+1))
        q109.angular.native.covariance.defect=defect;q109.radial.edge=lambda p,w:E(p+1)
        state.target.source.haar.rows=lambda n:dict(Ez=E(11))
        rows.current.current.profile.weighted_rho=lambda n:E(1)
        pm.operator=lambda n,p:(calls.append(('new_pressure',n,p)) or dict(KT=E(2),KA=E(3)))
        rows.bound=lambda n,b:dict(source_state=state,order=n,pole=pole+1,pre_fast_A1=E(17),fast_difference=dict(Ez=E(19)),slow_difference=dict(Ez=E(23)))
        bundle=dict(source_state=state,stage='post_first_sourced_wave_fixed_means',means_updated=False,order=12,covariance_pole=pole,change=[E(1)]*13)
        try:
            old=a.old(0,bundle)
            self.assertEqual([c for c in calls if c[0]=='old_signed'],[('old_signed',0),('old_signed',1),('old_signed',2)])
            self.assertIn(('old_pressure',1,pole),calls)
            expected=243*(old['initial']+(pole+1)*(q109.cv(old['signed_haar']['axial'])+E(12)))
            signature=lambda v:{k:a.terms for k,a in v.c.items()}
            self.assertEqual(signature(old['coefficient']),signature(expected))
            self.assertEqual(old['old_A1_haar'],0);self.assertTrue(old['old_A1_velocity_retained'])
            out=a.final(0,bundle)
            self.assertEqual(signature(out['coefficient']),signature(out['pre_fast']['pre_fast_haar']+(pole+2)*E(59)))
            self.assertIn(('new_pressure',1,pole),calls)
            for mutation in (dict(order=11),dict(means_updated=True),dict(source_state=object()),dict(stage='postslow')):
                with self.assertRaises(ValueError):a.final(0,dict(bundle,**mutation))
            q109.angular.native.covariance.defect=lambda st,n,**kw:dict(order=n,pole=pole+1,signed=E(1))
            with self.assertRaises(ValueError):a.old(0,bundle)
            q109.angular.native.covariance.defect=defect
            fresh.pressure=lambda n,p:dict(order=n,pole=p+1,KT=E(2),KA=E(3))
            with self.assertRaises(ValueError):a.old(0,bundle)
            a.state=object()
            with self.assertRaises(ValueError):a.final(0,bundle)
        finally:fresh.bound=oldbound;q109.angular.native.covariance.defect=olddefect;q109.radial.edge=oldedge
        with self.assertRaises(TypeError):q109.FullAxial(None)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
