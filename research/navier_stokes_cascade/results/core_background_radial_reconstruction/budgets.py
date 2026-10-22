"""Low-order full velocity/primitive norms from finalized F,U rows."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m


head=load('radial_reconstruction_head',HERE.parent/'core_background_finite_head_budget/budgets.py')
first=load('radial_reconstruction_first',HERE.parent/'global_first_coefficient_c3/global_first.py')
second=load('radial_reconstruction_second',HERE.parent/'global_second_coefficient_h23/global_second.py')
exact=head.exact


def positive_norm(x):
    x=exact(x)
    if x<0:raise ValueError('Nonnegative exact norm required')
    return x


def add(*polys):
    out={}
    for p in polys:
        for a,c in p.items():out[a]=out.get(a,Q(0))+c
    return {a:c for a,c in out.items() if c}


def scale(c,p):return {a:c*v for a,v in p.items() if c*v}


def mul(a,b,m):
    out={}
    for (r,k),c in a.items():
        for (s,l),d in b.items():
            if r+k+s+l<=m:out[r+s,k+l]=out.get((r+s,k+l),Q(0))+c*d
    return {a:c for a,c in out.items() if c}


def radius_series(m,lower,upper,inverse=False):
    out={(0,0):1/lower if inverse else upper};fall=Q(1)
    exponent=Q(-1,2) if inverse else Q(1,2)
    for j in range(1,m+1):
        fall*=exponent-j+1
        out[j,0]=abs(fall)*2**j*lower**((-1 if inverse else 1)-2*j)/factorial(j)
    return out


def reconstruct(h,n,m,Xhi,Rlo,Rhi,Frows,Urows):
    """Rows r bound |partial_X^r(F,U)|_{C_eta^(m+1-r)}.

    U rows hold on [0,Xhi], F rows on [Rlo^2/2,Xhi]. No primitive
    datum is reset at a velocity-zero interval. Only orders 0..2 used here.
    """
    h,Xhi,Rlo,Rhi=map(exact,(h,Xhi,Rlo,Rhi))
    if not 0<h<Q(1,100) or type(n) is not int or n<0 or type(m) is not int or not 0<=m<=2:
        raise ValueError('Retained h, nonnegative coefficient and low mixed order 0..2 required')
    if not 0<Rlo<=Rhi or not Rlo*Rlo<=2*Xhi<=Rhi*Rhi:
        raise ValueError('Nonempty positive annulus and valid square-root radius bounds required')
    Frows,Urows=tuple(map(positive_norm,Frows)),tuple(map(positive_norm,Urows))
    if len(Frows)!=m+1 or len(Urows)!=m+1:raise ValueError('Every radial row through m required')
    # Factorial angular rows convert to RAW bounds before bivariate Taylor use.
    uraw=lambda r,k:factorial(k)*Urows[r]
    U={(r,k):uraw(r,k)/Q(factorial(r)*factorial(k)) for r in range(m+1) for k in range(m+1-r)}
    F={(r,k):Frows[r]/Q(factorial(r)) for r,k in U}
    def primitive(r,k):return Xhi*uraw(0,k) if r==0 else uraw(r-1,k)
    M={(r,k):primitive(r,k)/Q(factorial(r)*factorial(k)) for r,k in U}
    Meta={(r,k):primitive(r,k+1)/Q(factorial(r)*factorial(k)) for r,k in U}
    eta={(0,0):Q(1),(0,1):Q(1)} if m else {(0,0):Q(1)}
    X={(0,0):Xhi,(1,0):Q(1)} if m else {(0,0):Xhi}
    d={(0,0):Q(1)}
    if m>=1:d[0,1]=Q(2)
    if m>=2:d[0,2]=Q(1)
    ell=1-2*h;Linv={(0,0):1/ell}
    for k in range(1,m+1):Linv[0,k]=(4*h*Linv[0,k-1]+(2*h*Linv[0,k-2] if k>=2 else 0))/ell
    times=lambda a,b:mul(a,b,m)
    D=Q(1,2)-h;lam=2*n*h
    V=times(Linv,add(scale(2,times(eta,times(X,U))),scale(2*(D+lam),times(eta,M)),times(d,Meta)))
    R=radius_series(m,Rlo,Rhi);Ri=radius_series(m,Rlo,Rhi,True)
    polys=dict(E=times(R,F),U=U,W=times(Ri,V),J=scale(2,times(Ri,times(Linv,times(eta,M)))))
    raw={name:{a:poly.get(a,Q(0))*factorial(a[0])*factorial(a[1]) for a in U} for name,poly in polys.items()}
    return dict(order=m,coefficient=n,raw=raw,
                through={j:max((v for rows in raw.values() for a,v in rows.items() if sum(a)<=j),default=Q(0)) for j in range(m+1)},
                primitive_raw={a:primitive(*a) for a in U},
                angular_input_order=m+1,
                scope='conditional full-row reconstruction on the same annulus; full incoming primitive and pressure-finalized repair retained')


FIRST_KEYS={'nu','x_velocity','leading_u5','x_cut','local_seed','leading_f3','leading_u3','rho','lam','geometric_k','e_star'}


def finalized_first(h,inner,cutoff,first_data,p_lower,rp_lower):
    if set(first_data)!=FIRST_KEYS:raise ValueError('Complete original first source/moment/linear repair data required')
    first_data={k:exact(v) for k,v in first_data.items()}
    inner,cutoff=tuple(map(positive_norm,inner)),tuple(map(positive_norm,cutoff))
    if len(cutoff)!=4 or cutoff[0]!=1:raise ValueError('Same eta-independent coefficient cutoff with unit value bound required')
    if len(inner)!=4 or exact(inner[0])!=exact(first_data['local_seed']):raise ValueError('Same whole-inner Ceta3 value seed and four radial rows required')
    repair=first.finalized_k13(h=h,**first_data)
    jets=second.finalized_first_jets(inner,cutoff,repair['repair']['coefficient_Ceta3'],p_lower,rp_lower)
    return dict(repair=repair,jets=jets)


def low_comparison(h,Xhi,Rlo,Rhi,leading_F,leading_U,first_state,cutoff_Euler2,geometry):
    leading=reconstruct(h,0,2,Xhi,Rlo,Rhi,leading_F,leading_U)
    current=reconstruct(h,1,1,Xhi,Rlo,Rhi,first_state['jets']['F'][:2],first_state['jets']['U'][:2])
    K=positive_norm(cutoff_Euler2)
    if K<1:raise ValueError('Actual Borel cutoff Euler norm through two, padded by one, required')
    def lead(m):
        if m not in (1,2):raise ValueError('Only the proved low-order leading inputs are supplied')
        return leading['through'][m]
    def profile(n,m):
        if (n,m)!=(1,1):raise ValueError('Only the proved coefficient-one first-jet input is supplied')
        return current['through'][1]
    def cutoff(m):
        if m!=2:raise ValueError('Only the required Borel cutoff order two is supplied')
        return K
    # Rlo,Rhi here are NATIVE radii. Band radii in geometry are separate.
    g=dict(geometry);band_lo=g.pop('Rmin');band_hi=g.pop('Rmax')
    band_lo,band_hi=map(exact,(band_lo,band_hi))
    if not 0<band_lo<=band_hi or 4*exact(Rlo)**2>band_lo**2 or exact(Xhi)<2*band_hi**2:
        raise ValueError('Native annular bounds must cover the ENTIRE band chart for s in [1/4,4]')
    b=head.BackgroundBudget(h,lead,profile,cutoff,band_lo,band_hi,'euler')
    return dict(leading=leading,first=current,comparison=b.comparison(**g),
                scope='explicit low-order Cphi input join; normalized Borel schedule and actual selector value substitution remain premises')


def fixture():
    h=Q(1,200);inner=(Q(2),Q(3),Q(5),Q(7))
    data=dict(nu=Q(1),x_velocity=Q(3),leading_u5=Q(5),x_cut=Q(1),local_seed=inner[0],
              leading_f3=Q(7),leading_u3=Q(5),rho=Q(2),lam=Q(1,100),geometric_k=Q(1),e_star=Q(1))
    cutoff=second.coefficient_cutoff_jets(Q(1,8),Q(1,2))
    state=finalized_first(h,inner,cutoff,data,Q(1,2),Q(1))
    geometry=dict(Rmin=Q(1),Rmax=Q(2),U=Q(2),lambda_max=Q(3),gmin=Q(1,2),Llow=Q(1),Lhigh=Q(2),Dbox=Q(1))
    result=low_comparison(h,Q(8),Q(1,2),Q(4),(Q(7),Q(11),Q(13)),(Q(5),Q(17),Q(19)),state,Q(3),geometry)
    return dict(h=h,inner=inner,first_data=data,cutoff=cutoff,state=state,result=result,geometry=geometry)
