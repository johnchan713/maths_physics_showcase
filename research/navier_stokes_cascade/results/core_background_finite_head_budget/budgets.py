"""Finite-head background constants and an explicit full-path phase comparison."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module);return module


phase=load("finite_head_phase",HERE.parent/"core_phase_coefficient_budget/budgets.py")
P,order,exact=phase.P,phase.order,phase.exact


def nonnegative(value):
    value=exact(value)
    if value<0:raise ValueError("Nonnegative finite bound required")
    return value


def increment_constant(h,n,m,profile,cutoff_euler):
    h=exact(h);order(m)
    if not 0<h<Q(1,100) or type(n) is not int or n<1:raise ValueError("Fixed construction h and positive coefficient order required")
    profile,cutoff=nonnegative(profile),nonnegative(cutoff_euler)
    if cutoff<1:raise ValueError("Padded Euler cutoff norm through m+1 required")
    return max(Q(1),2*cutoff*(1+2*n*h)**m*profile)


def finite_head(m,head):
    order(m);N=max(2,m)
    if set(head)!=set(range(1,N)):raise ValueError("Every finite-head coefficient and no tail values required")
    values={n:nonnegative(c) for n,c in head.items()}
    tail=Q(1,2**(N-1))
    return dict(start=N,head=sum(values.values(),Q(0)),tail=tail,correction=sum(values.values(),Q(0))+tail)


def chart_constants(m,Rmin,Rmax):
    order(m);lo,hi=exact(Rmin),exact(Rmax)
    if lo<=0 or lo>hi:raise ValueError("Positive closed annular radius bounds required")
    Y=sum((Q(32**i*64**j,3) for i in range(m+1) for j in range(m+1-i) if i+j),Q(0))
    exp=lambda a:sum((a**k/Q(factorial(k)) for k in range(m+1)),Q(0))
    E,Eh=exp(Y),exp(Y/2)
    # Zero-constant increments of X and eta; s in [1/4,4].
    Rpoly=hi*hi+(2*hi if m>=1 else 0)+(1 if m>=2 else 0)
    X=2*(Rpoly*E-hi*hi)
    eta=(1+(4 if m>=1 else 0))*Eh-1
    C=sum((Y**i*X**a*eta**b/Q(factorial(i)*factorial(a)*factorial(b))
           for i in range(m+1) for a in range(m+1-i) for b in range(m+1-i-a)),Q(0))
    reciprocal=sum((lo**(-j-1) for j in range(m+1)),Q(0))
    return dict(log=Y,delta_X=X,delta_eta=eta,composition=C,
                F=factorial(m)*4*E*reciprocal*C,
                G=factorial(m)*4*E*C,btilde=factorial(m)*2*Eh*C)


class BackgroundBudget:
    """No cutoff schedule is selected or changed by this finite norm recipe."""
    def __init__(self,h,leading,profiles,cutoff,Rmin,Rmax,derivative_basis,trace=None):
        self.h=exact(h)
        if not 0<self.h<Q(1,100):raise ValueError("Fixed construction h required")
        if derivative_basis!="euler":raise ValueError("Use Euler q derivatives, not frozen q^j partial_q^j with identical constants")
        self.leading,self.profiles,self.cutoff=leading,profiles,cutoff
        self.Rmin,self.Rmax=exact(Rmin),exact(Rmax);chart_constants(0,self.Rmin,self.Rmax)
        self.trace=trace

    def leading_bound(self,m):
        order(m)
        if self.trace is not None:self.trace.add(("leading",m))
        chart=chart_constants(m,self.Rmin,self.Rmax)
        return (chart["F"]+chart["G"])*nonnegative(self.leading(m))

    @lru_cache(None)
    def bound(self,m):
        order(m)
        if self.trace is not None:self.trace.add(("leading",m));self.trace.add(("cutoff",m+1))
        L=nonnegative(self.leading(m));K=nonnegative(self.cutoff(m+1))
        head={}
        for n in range(1,max(2,m)):
            if self.trace is not None:self.trace.add(("profile",n,m))
            head[n]=increment_constant(self.h,n,m,self.profiles(n,m),K)
        out=finite_head(m,head);chart=chart_constants(m,self.Rmin,self.Rmax)
        total=L+out["correction"];fg=chart["F"]+chart["G"]
        return dict(F=chart["F"]*total,G=chart["G"]*total,btilde=chart["btilde"]*total,
                    leading_FG=fg*L,full_FG=fg*total,error_FG=4*fg*out["correction"],
                    head=out,chart=chart,cutoff_input_order=m+1)

    def phase_inputs(self,path,Lmin,U,lambda_max,gmin):
        cb=phase.carrier_bounds(lambda_max,U,gmin,self.Rmax,Lmin)
        return phase.PhaseBudget(F=lambda j:P(self.bound(j)["F"]),G=lambda j:P(self.bound(j)["G"]),
            btilde=lambda j:P(self.bound(j)["btilde"]),path=path,Rmin=self.Rmin,Rmax=self.Rmax,
            Lmin=Lmin,u=U,lambda_max=lambda_max,carriers=cb)

    def comparison(self,U,lambda_max,gmin,Llow,Lhigh,Dbox):
        U,lam,g,lo,hi,D=map(exact,(U,lambda_max,gmin,Llow,Lhigh,Dbox))
        if min(lam,g,lo,hi)<=0 or min(U,D)<0 or hi<lo:raise ValueError("Valid leading and whole-path geometry constants required")
        B=max(Q(1),lam);cb=phase.carrier_bounds(lam,U,g,self.Rmax,lo)
        p,pz=cb["theta"],cb["axial"];p0=p-1
        lead1=self.leading_bound(1);lead2=self.leading_bound(2)
        err1=self.bound(1)["error_FG"];full1=self.bound(1)["full_FG"]
        groups=dict(tilt=2*B*U/(lo*g),radial_box=hi*(p+pz)*3*D*lead2,
                    radial_background=hi*(p+pz)*err1,radial_rounding=hi*lead1,
                    angular_rounding=1/self.Rmin,angular_box=p0*D/self.Rmin**2,
                    axial_transport=hi*(p+pz)*full1)
        return dict(Cphi=sum(groups.values(),Q(0)),groups=groups,carriers=cb,
                    retained_dependencies="leading jets through two; finalized coefficient-one tuple jets through one; Borel cutoff Euler jets through two; normalized diagonal schedule")


def length_compatibility(lambda_min,U,Cphi,radius,Lreq,slow_cost,slow_tolerance):
    lm,U,C,r,L,cs,tau=map(exact,(lambda_min,U,Cphi,radius,Lreq,slow_cost,slow_tolerance))
    if min(lm,r,L,cs,tau)<=0 or min(U,C)<0 or cs<1 or tau>1:
        raise ValueError("Retained positive selector data and its original slow ranges required")
    beta=lm/(4*(1+U*U)**2)
    return dict(beta=beta,inverse_normal_square=4/beta,
                gates=dict(original_slow_selector_implies_rounding_gate=tau/cs<=1,
                           original_length_budget_dominates_normal_margin=beta*L*L>=144*r*r*C*C),
                scope="conditional on all original selector checks and actual constants; do not reselect Lreq or a physical band")


def fixture(trace=None):
    return BackgroundBudget(Q(1,200),leading=lambda m:Q(m+2),profiles=lambda n,m:Q(n+m+1),
                            cutoff=lambda m:Q(m+1),Rmin=1,Rmax=2,derivative_basis="euler",trace=trace)
