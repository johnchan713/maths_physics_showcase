import unittest,json
from fractions import Fraction
import sympy as s
import audit,q106

class Tests(unittest.TestCase):
    def test_full_rows_average_to_ledger_with_correlations_and_viscosity(self):
        R,Z,T,y=s.symbols('R Z T y',real=True);e=s.Rational(1,3);h=(R-1)**3*(2-R)**3
        avg=lambda f:s.integrate(s.expand(f),(y,-s.pi,s.pi))/(2*s.pi)
        Dr=lambda f:s.diff(f,R)+2*R*s.diff(f,y)
        De=lambda f,k:Dr(f)+k*f/R
        Dz=lambda f:e*s.diff(f,Z)
        ts=lambda f:-e*s.diff(f,T)+2*s.diff(f,y)
        x=h*(Z+s.cos(y));u=h*(T+s.cos(y)+s.sin(y));z=h*(1+Z**2+s.sin(y))
        base=dict(b=e*R,V=1+R,G=R*Z);b,V,G=(base[k] for k in ('b','V','G'))
        W={k:h*(i+2+s.cos(y)) for i,k in enumerate(('rr','rt','rz','tt','zt','zz'))}
        st=dict(t=h*s.cos(y),z=h*(Z+s.sin(y)));pm=h*(R*Z+s.sin(y))
        full_t=ts(u)+De(b*u+V*x+x*u+W['rt']-st['t'],2)+Dz(G*u+V*z+z*u+W['zt'])-e*(Dr(Dr(u))+Dr(u)/R-u/R**2+Dz(Dz(u)))
        full_z=ts(z)+De(b*z+G*x+x*z+W['rz']-st['z'],1)+Dz(2*G*z+z*z+W['zz']+pm)-e*(De(Dr(z),1)+Dz(Dz(z)))
        means={k:avg(f) for k,f in dict(x=x,u=u,z=z).items()}
        products={k:avg(f) for k,f in dict(xx=x*x,uu=u*u,zz=z*z,xu=x*u,zu=z*u,xz=x*z).items()}
        out=q106.rows(R,Z,T,e,means,products,{k:avg(f) for k,f in W.items()},base,avg(pm),{k:avg(f) for k,f in st.items()})
        self.assertEqual(s.simplify(avg(full_t)-out['theta']),0)
        self.assertEqual(s.simplify(avg(full_z)-out['axial']),0)
        wrong=dict(products,xu=means['x']*means['u'])
        bad=q106.rows(R,Z,T,e,means,wrong,{k:avg(f) for k,f in W.items()},base,avg(pm),{k:avg(f) for k,f in st.items()})
        self.assertNotEqual(s.simplify(out['theta']-bad['theta']),0)
    def test_positive_covariance_all_moments_zero_but_pointwise_row_nonzero(self):
        R,Z,T=s.symbols('R Z T',real=True);h=(R-1)**3*(2-R)**3;f=h*h
        zero=dict(x=0,u=0,z=0);prod={k:0 for k in ('xx','uu','zz','xu','zu','xz')}
        W=dict(rr=f,tt=f,zz=f,rt=f,rz=0,zt=0);base=dict(b=0,V=0,G=0)
        split=q106.radial_split(R,Z,T,1,zero,prod,W,base)
        self.assertEqual(split['r0'],0);self.assertEqual(split['K'],-f)
        out=q106.rows(R,Z,T,1,zero,prod,W,base,-f,dict(t=0,z=0))
        self.assertEqual(s.simplify(out['theta']-s.diff(f,R)-2*f/R),0)
        self.assertEqual(out['axial'],0);self.assertNotEqual(out['theta'],0)
        # Exact Q102 covariance moment formulas, including axial pressure part.
        P=s.integrate((W['tt']-W['rr'])/R,(R,1,2))
        Jt=s.integrate(R*W['zt'],(R,1,2))
        Jz=s.integrate(R*W['zz']-R*(W['rr']+W['tt'])/2,(R,1,2))
        self.assertEqual((P,Jt,Jz),(0,0,0))
        self.assertEqual(s.integrate(s.expand(R**2*out['theta']),(R,1,2)),0)
        M=s.Matrix([[1,1,0],[1,1,0],[0,0,1]])
        self.assertEqual(sorted(M.eigenvals()),[0,1,2])
        y=s.symbols('y',real=True);v=s.sqrt(2)*h*s.Matrix([s.cos(y),s.cos(y),s.sin(y)])
        avg=lambda f:s.integrate(s.expand_trig(f),(y,-s.pi,s.pi))/(2*s.pi)
        self.assertEqual((v*v.T).applyfunc(avg),f*M)
    def test_moving_pressure_haar_commutes_but_memory_survives(self):
        R,r,y,Z=s.symbols('R r y Z',real=True);h=(R-1)**3*(2-R)**3;t=R-1
        chi=3*t*t-2*t**3+Z*t*t*(1-t)**2
        f=h*(s.cos(y-2*R)+Z);K=h*(1+s.sin(y-2*R))
        def integ(g,hi):return s.integrate(s.expand(g.subs({R:r,y:y+2*(r-R)},simultaneous=True)),(r,1,hi))
        T0=integ(f,R)-chi*integ(f,2)
        avg=lambda g:s.integrate(s.expand(g),(y,-s.pi,s.pi))/(2*s.pi)
        expected=integ(h*Z,R)-chi*integ(h*Z,2)
        self.assertEqual(s.simplify(avg(K+T0)-h-expected),0)
        self.assertNotEqual(s.simplify(avg(K+T0)-expected),0)
        self.assertNotEqual(s.simplify(s.diff(expected,Z)-(integ(h,R)-chi*integ(h,2))),0)
    def test_radial_split_keeps_averaged_products(self):
        R,Z,T=s.symbols('R Z T',positive=True);e=s.Rational(1,3)
        m=dict(x=R*Z,u=R*T,z=R*Z**2);p=dict(xx=R**2*Z**2+1,uu=R**2*T**2+2,xz=R**2*Z**3+3)
        W=dict(rr=R**2,rz=R*Z,tt=2*R**2)
        a=q106.radial_split(R,Z,T,e,m,p,W,dict(b=R,V=2*R,G=R*Z))
        direct=e*s.diff(m['x'],T)-s.diff(a['S'],R)-a['S']/R-e*s.diff(a['H'],Z)+a['C']/R+e*(s.diff(m['x'],R,2)+s.diff(m['x'],R)/R-m['x']/R**2+e**2*s.diff(m['x'],Z,2))
        self.assertEqual(s.simplify(direct-a['r0']-s.diff(a['K'],R)),0)
    def test_full_join_no_extra_projection_factor(self):
        self.assertEqual(q106.join(dict(theta=2,axial=3),dict(Etheta=5,Ez=7),dict(Etheta=11,Ez=13),17),dict(theta=18,axial=40))
        self.assertFalse(q106.record()['actual_fine_haar_norm_supplied'])
    def test_positive_envelope_and_domain_guards(self):
        Q=Fraction
        m=dict(u=[2,3,5],z=[7,11,13]);f=dict(At=[17,19],Bt=[23,29],Az=[31,37],Bz=[41,43])
        a=q106.envelope(0,Q(1,3),m,f,2,4)
        self.assertEqual(a['theta'],Q(1,3)*(3+10+6+8+Q(10,9))+19+4*17+Q(29,3))
        self.assertEqual(a['axial'],Q(1,3)*(11+26+22+Q(26,9))+37+2*31+Q(43,3))
        for n,e in ((-1,1),(0,0),(0,2),(True,1)):
            with self.assertRaises(ValueError):q106.envelope(n,e,m,f,2,4)
        with self.assertRaises(ValueError):q106.envelope(1,1,m,f,2,4)
        with self.assertRaises(ValueError):q106.envelope(0,1,m,dict(f,At=[-1,2]),2,4)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
