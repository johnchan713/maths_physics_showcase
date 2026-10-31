"""Exact radial adjoints in the native volume/shift convention."""
import sympy as s


def weights(r,c,V,G):
    w=c-r**2/2
    a=2*V*w/r
    return dict(zz=r,rr=s.diff(w,r)-w/r,tt=w/r,
                rz=2*r*s.diff(G,r),rt=s.diff(a,r)-2*a/r,
                radial_viscosity=s.simplify(s.diff(w,r,2)-s.diff(w,r)/r))


def record():
    r,c=s.symbols('R c',positive=True)
    d=weights(r,c,s.Function('V')(r),s.Function('G')(r))
    assert s.simplify(d['zz']+d['tt']+d['rr'])==0
    assert d['radial_viscosity']==0
    from fractions import Fraction as F
    k=F(1,100000)
    return dict(components=list(d),isotropic_response=0,radial_viscous_response=0,
                remainder=str(1+min(2-2*k,2-k,2)),target=str(3-5*k))
