"""Exact Haar mean-row ledger; no fine native norm claim."""
import sympy as s
from fractions import Fraction as Q

def envelope(n,epsilon,mean,flux,r,r2):
    """Conditional inequality (6), not native coefficient certification."""
    if type(n) is not int or n<0:raise ValueError('Nonnegative integer order required')
    e=Q(epsilon);r=Q(r);r2=Q(r2)
    if not 0<e<=1 or min(r,r2)<0:raise ValueError('Positive epsilon and nonnegative norm coefficients required')
    if set(mean)!={'u','z'} or set(flux)!={'At','Bt','Az','Bz'}:raise ValueError('Complete Haar flux channels required')
    def get(a,j):
        if len(a)<=j:raise ValueError('Insufficient Haar coefficient jets')
        v=Q(a[j])
        if v<0:raise ValueError('Nonnegative norm coefficients required')
        return v
    d=n+1;dd=d*(d+1)
    def field(k,radial):
        a=mean[k];a0=get(a,n);a1=get(a,n+1);a2=get(a,n+2)
        return e*(d*a1+dd*a2+r*d*a1+radial*r2*a0+e*e*dd*a2)
    at,bt,az,bz=(flux[k] for k in ('At','Bt','Az','Bz'))
    return dict(theta=field('u',1)+d*get(at,n+1)+2*r*get(at,n)+e*d*get(bt,n+1),
                axial=field('z',0)+d*get(az,n+1)+r*get(az,n)+e*d*get(bz,n+1))

def rows(R,Z,T,eps,means,products,wave,base,pressure,stress):
    """All arguments are already Haar averaged, except fixed base factors.

    products MUST be averages of products, not products of averaged fields.
    pressure is the complete averaged shifted pressure, including local memory.
    Label parameters and base multipliers are auxiliary-independent.
    """
    x,u,z=(means[k] for k in ('x','u','z'))
    b,V,G=(base[k] for k in ('b','V','G'))
    D=lambda f:s.diff(f,R)
    Dz=lambda f:eps*s.diff(f,Z)
    De=lambda f,e:D(f)+e*f/R
    rt=b*u+V*x+products['xu']+wave['rt']-stress['t']
    zt=G*u+V*z+products['zu']+wave['zt']
    rz=b*z+G*x+products['xz']+wave['rz']-stress['z']
    zz=2*G*z+products['zz']+wave['zz']+pressure
    vt=D(D(u))+D(u)/R-u/R**2+Dz(Dz(u))
    vz=De(D(z),1)+Dz(Dz(z))
    return dict(theta=-eps*s.diff(u,T)+De(rt,2)+Dz(zt)-eps*vt,
                axial=-eps*s.diff(z,T)+De(rz,1)+Dz(zz)-eps*vz)

def radial_split(R,Z,T,eps,means,products,wave,base):
    x,u,z=(means[k] for k in ('x','u','z'))
    b,V,G=(base[k] for k in ('b','V','G'))
    S=2*b*x+products['xx']+wave['rr']
    H=b*z+G*x+products['xz']+wave['rz']
    C=2*V*u+products['uu']+wave['tt']
    r0=eps*s.diff(x,T)-S/R-eps*s.diff(H,Z)+C/R+eps**3*s.diff(x,Z,2)
    K=-S+eps*(s.diff(x,R)+x/R)
    return dict(S=S,H=H,C=C,r0=r0,K=K)

def join(pre,fast,slow,tail):
    """Positive coefficients, same norm/pole; full differences need no Q factor."""
    return dict(theta=pre['theta']+fast['Etheta']+slow['Etheta'],
                axial=pre['axial']+tail+fast['Ez']+slow['Ez'])

def record():
    return dict(status='exact-pointwise-haar-ledger-and-moment-blind-covariance-control',
                products_averaged_before_factorization=True,
                full_local_pressure_and_both_axial_viscosities_retained=True,
                all_three_global_moments_can_vanish_in_comparison_covariance=True,
                canonical_wave_counterexample=False,actual_fine_haar_norm_supplied=False,
                full_cone_residual_pulse_or_blowup_certified=False)
