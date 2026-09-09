"""Couple scheduled pressure to B.15 using a second coefficient assembly.

The earlier pilot expands R1/R2. Here the original W, H, U, Pi products are
expanded directly. Agreement on the old seed is a useful implementation check.
The old source files and their frozen evidence are not modified.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
import mpmath as mp

PREVIOUS = Path(__file__).resolve().parent.parent/'nonlinear_axis_pilot'
sys.path.insert(0,str(PREVIOUS))
from inner import (Parameters,InnerProfile,cut,add,scale,mul,derivative,
                   product_n,axis_series as seed_axis_series,inner_moments)


@dataclass(frozen=True)
class CoupledParameters(Parameters):
    datum: object

    @classmethod
    def from_schedule(cls,schedule,data):
        return cls(schedule.h,mp.mpf(data['j']),mp.mpf(data['sigma']),
                   mp.mpf(data['Lambda_over_P_squared'])*mp.exp(2*schedule.log_p),
                   mp.mpf(0),schedule)


class SeedDatum:
    """Only for the independent regression against the previous R1/R2 solver."""
    def __init__(self,amplitude):
        self.amplitude=amplitude

    def pressure_jet(self,eta,degree):
        from schedule import power_jet
        return [-self.amplitude*v for v in power_jet(eta,2,degree)]

    def pressure(self,eta):
        return -self.amplitude/(1+eta*eta)**2


def axis_series(p,eta,degree):
    base=seed_axis_series(p,eta,degree)
    pressure=p.datum.pressure_jet(eta,degree+1)
    # p.pressure=0 in the seed helper, so replace only its pressure contribution.
    correction=add(scale(mul(base['d'],derivative(pressure),degree),-1,degree),
                   scale(mul(base['eta'],pressure,degree),4*p.A,degree),m=degree)
    base['pressure']=pressure
    base['Z']=add(base['Z'],correction,m=degree)
    return base


class CoupledProfile(InnerProfile):
    def point(self,Y,eta=None):
        eta=self.eta0 if eta is None else eta
        p=self.parameters
        logg=self.log_g0 if eta==self.eta0 else p.log_g(eta)
        return (self.value(self.phi,Y,eta),4*eta+p.j+self.value(self.u,Y,eta)/p.lam,
                p.datum.pressure(eta)+self.value(self.pressure_increment,Y,eta)/p.lam,mp.exp(logg))


def construct(p,eta0,degree,extra_eta=3):
    if not isinstance(degree,int) or not 2<=degree<=40 or not 2<=extra_eta<=6 or abs(eta0)>1:
        raise ValueError('Invalid degree, derivative budget or eta center')
    total=degree+extra_eta
    base=axis_series(p,eta0,total)
    phi=[cut([mp.mpf(1)],total)]
    u=[cut([mp.mpf(0)],total)]
    pressure=[cut([mp.mpf(0)],total)]
    logg=p.log_g(eta0)
    g2=[mp.exp(2*logg)]
    for k in range(total):
        g2.append(2*p.lam*mp.fsum(base['zeta'][i]*g2[k-i] for i in range(k+1))/(k+1))
    for n in range(degree):
        m=total-n-1
        eta,d,invL=[base[k] for k in ('eta','d','inverse_L')]
        # Actual U and Pi coefficients, with one extra eta derivative available.
        U=[scale(row,1/p.lam,m+1) for row in u]
        U[0]=add(U[0],base['U'],m=m+1)
        Pi=[scale(row,1/p.lam,m+1) for row in pressure]
        Pi[0]=add(Pi[0],base['pressure'],m=m+1)
        H=[mul(d,row,m) for row in U]
        H[0]=add(H[0],scale(eta,p.D,m),m=m)
        W=[scale(add(scale(mul(eta,row,m),2*p.D,m),mul(d,derivative(row),m),m=m),
                 -1/mp.mpf(i+1),m) for i,row in enumerate(U)]
        W[0]=add(W[0],[mp.mpf(1)],m=m)
        Hphi=product_n(H,phi,n,m)
        angular=add(product_n(W,phi,n,m),product_n(W,phi,n,m,right_radial_weight=True),
                    scale(phi[n],p.h,m),scale(mul(eta,product_n(U,phi,n,m),m),-2*p.h,m),
                    scale(mul(base['zeta'],Hphi,m),p.lam,m),
                    product_n(H,phi,n,m,right_derivative=True),m=m)
        axial=add(scale(product_n(W,u,n,m,right_radial_weight=True),1/p.lam,m),
                  scale(U[n],p.A,m),scale(mul(eta,product_n(U,U,n,m),m),-2*p.A,m),
                  product_n(H,U,n,m,right_derivative=True),mul(d,derivative(Pi[n]),m),
                  scale(mul(eta,Pi[n],m),-4*p.A-2*n,m),m=m)
        phi.append(scale(mul(invL,angular,m),1/(2*p.lam*(n+1)*(n+2)),m))
        u.append(scale(mul(invL,axial,m),1/mp.mpf(2*(n+1)**2),m))
        pressure.append(scale(mul(g2,product_n(phi,phi,n,m),m),1/mp.mpf(n+1),m))
    return CoupledProfile(p,eta0,degree,phi,u,pressure,logg)


def diagnostics(profile,Y):
    """Evaluate the scalar source identities at a point after construction."""
    p,e=profile.parameters,profile.eta0
    d,L=1-e*e,1-2*p.h*e*e
    Phi,U,Pi,g=profile.point(Y)
    val=profile.value
    py,pyy,pe=[val(profile.phi,Y,**kw) for kw in
               ({'radial_order':1},{'radial_order':2},{'eta_order':1})]
    uy,uyy,ue=[val(profile.u,Y,**kw) for kw in
               ({'radial_order':1},{'radial_order':2},{'eta_order':1})]
    average=[scale(row,1/mp.mpf(n+1),len(row)-1) for n,row in enumerate(profile.u)]
    avgU=4*e+p.j+val(average,Y)/p.lam
    avgUe=4+val(average,Y,eta_order=1)/p.lam
    W=1-2*p.D*e*avgU-d*avgUe
    H=p.D*e+d*U
    Pi_e=p.datum.pressure_jet(e,1)[1]+val(profile.pressure_increment,Y,eta_order=1)/p.lam
    PiX=val(profile.pressure_increment,Y,radial_order=1)
    angular_terms=[-W*(1+Y*py/Phi),-p.h*(1-2*e*U),-H*(p.lam*p.zeta(e)+pe/Phi)]
    axial_terms=[-W*Y*uy/p.lam,-p.A*(1-2*e*U)*U,-H*(4+ue/p.lam),
                 -d*Pi_e,4*p.A*e*Pi,2*e*Y*PiX/p.lam]
    angular_lhs=-2*L*p.lam*(Y*pyy+2*py)/Phi
    axial_lhs=-2*L*(Y*uyy+uy)
    def balance(lhs,terms):
        denom=max(mp.mpf(1),abs(lhs),*(abs(t) for t in terms))
        return abs(lhs-mp.fsum(terms))/denom
    return dict(Y=Y,Phi=Phi,U=U,Pi=Pi,log_g=profile.log_g0,
                angular_residual=balance(angular_lhs,angular_terms),
                axial_residual=balance(axial_lhs,axial_terms),
                pressure_residual=abs(PiX-(g*Phi)**2)/max(abs(PiX),abs((g*Phi)**2)),
                normalized_axial_pressure_contribution=abs(-d*Pi_e+4*p.A*e*Pi)/max(mp.mpf(1),abs(axial_lhs),*(abs(t) for t in axial_terms)),
                angular_slope=mp.mpf('.5')+Y*py/Phi,
                moments=inner_moments(profile,Y))


def compare(a,b,ys):
    errors=[]
    for y in ys:
        for name in ('phi','u'):
            for dy,de in ((0,0),(1,0),(2,0),(0,1),(0,2),(0,3)):
                x=a.value(getattr(a,name),y,radial_order=dy,eta_order=de)
                z=b.value(getattr(b,name),y,radial_order=dy,eta_order=de)
                errors.append(abs(x-z)/max(mp.mpf(1),abs(x),abs(z)))
    return max(errors)


def full_residual(profile,point):
    """Reuse Cartesian derivative arithmetic and supply the new pressure gradient.

    The old field evaluator sees p.pressure=0 and therefore includes only the
    pressure increment. The missing axis datum contributes only to p_z.
    """
    from physical import field,q_value
    p=profile.parameters
    u,pressure=field(profile,point)
    q=q_value(point[2],point[3],p)
    eta=point[2]*q**(-p.D)
    pi0,pie=p.datum.pressure_jet(eta,1)
    pressure.g[2]+=q**(-2*p.A-p.D)*((1-eta*eta)*pie-4*p.A*eta*pi0)/(1-2*p.h*eta*eta)
    terms=[[v.g[3],mp.fsum(u[j].v*v.g[j] for j in range(3)),
            -mp.fsum(v.H[j][j] for j in range(3)),pressure.g[i]] for i,v in enumerate(u)]
    totals=[mp.fsum(row) for row in terms]
    scales=[max(abs(x) for x in row) for row in terms]
    leading=[totals[i]+u[i].H[2][2] for i in (1,2)]
    return dict(terms=terms,total=totals,scales=scales,
                normalized_full=[abs(v)/s for v,s in zip(totals,scales)],
                normalized_leading=[abs(v)/scales[i] for v,i in zip(leading,(1,2))],
                divergence=abs(mp.fsum(u[i].g[i] for i in range(3)))/max(mp.mpf(1),*(abs(u[i].g[i]) for i in range(3))))


def physical_finite_difference(profile,point,step):
    """Old component-separated stencil plus an independently differenced Pi0.

    Only p_z changes. This deliberately shares the coordinate map and finite
    profile; it is a derivative cross-check, not an independent whole solver.
    """
    from physical import finite_difference,q_value
    out=finite_difference(profile,point,step)
    p=profile.parameters
    q=q_value(point[2],point[3],p)
    stiffness=max(mp.mpf(1),abs(p.lam*p.zeta(profile.eta0)),mp.sqrt(abs(p.lam*mp.diff(p.zeta,profile.eta0))))
    dz=step*q**p.D/stiffness
    values=[]
    for k in (-2,-1,1,2):
        z=point[2]+k*dz
        qq=q_value(z,point[3],p)
        values.append(qq**(-2*p.A)*p.datum.pressure(z*qq**(-p.D)))
    a,b,c,d=values
    out[2]+=(a-8*b+8*c-d)/(12*dz)
    return out
