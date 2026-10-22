"""Conservative complete mean-size closure and the first actual source recipe."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


haar=load("source_size_haar",HERE.parent/"core_initial_haar_row_budget/budgets.py")
P=haar.P;order=haar.order;exact=haar.exact
STATE=("x","u","z")
FIXED=("b","V","G","sigma","graph","c","r","rho","Kp","Kmom")
LINEAR=("Eik","bR","Bz","normal","K","cutoff")


def pole_ledger(p):
    """Initial inputs at p; slow response has closed-interior weighted bounds."""
    p=exact(p)
    if p<0:raise ValueError("Nonnegative initial finite pole required")
    return dict(initial=p,covariance=2*p,pre_fast_rows=2*p,post_fast_means=2*p,
                fresh_post_fast_rows_and_moments=4*p,slow_response=Q(0),final_means=2*p,
                wave_wave=2*p,mean_wave=3*p,source=3*p)


def require(data,keys):
    if set(data)!=set(keys):raise ValueError("All and only the declared full-state inputs required")
    if all(isinstance(v,P) for v in data.values()):return data
    if all(type(v) in (int,Q) for v in data.values()):
        if min(data.values())<0:raise ValueError("Nonnegative norm bounds required")
        return {k:exact(v) for k,v in data.items()}
    raise TypeError("Use all positive polynomials or all exact rational bounds")


def full_mean(n,data):
    """Through n full graph rows, pressure and scalars; field inputs at n+3."""
    order(n);d=Q(n+4)
    data=require(data,STATE+FIXED+("w",))
    x,u,z,b,V,G,sigma,g,c,r,rho,kp,km,w=(data[k] for k in STATE+FIXED+("w",))
    dr=g*d;time=d*(1+c);vis=d*d*(g*g+1)+r*dr+r*r
    srr=2*b*x+x*x+w;srz=b*z+G*x+x*z+w;stt=2*V*u+u*u+w
    gr=time*x+(dr+r)*srr+d*srz+r*stt+vis*x
    defect=km*gr;pm=kp*(gr+rho*defect)
    rt=b*u+V*x+x*u+w+sigma;rz=b*z+G*x+x*z+w+sigma
    at=G*u+V*z+z*u+w;az=2*G*z+z*z+w+pm
    et=time*u+(dr+2*r)*rt+d*at+vis*u
    ez=time*z+(dr+r)*rz+d*az+vis*z
    jt=km*(G*u+V*z+z*u+w)
    jz=km*(2*G*z+z*z+w+gr)
    return dict(gr=gr,P=defect,pm=pm,Etheta=et,Ez=ez,Jtheta=jt,Jz=jz,
                field_order=n+3,pressure_source_order=n+1)


def source_recipe(n,data):
    """Factorial through-order norms; coefficient inputs through n+3."""
    order(n);d=Q(n+4)
    keys=("b","G","graph","c","r")+LINEAR+("T","R","Pi","U","M")
    data=require(data,keys)
    b,G,g,c,r,e,bR,Bz,N,K,psi,T,R,pi,U,M=(data[k] for k in keys)
    dr=g*d
    slow=d+b*dr+G*d
    phase=2*e
    base=bR+Bz+b*r
    viscosity=d*d*(g*g+1)+r*dr+r*r
    mixed=4*N*(dr+d)+2*N*(dr+d+r)
    connection=4*N*r
    E=slow+phase+base+viscosity+mixed+connection
    pressure=(dr+d)*pi
    # Principal homogeneous cancellation applies to T, not to the curl R.
    curl=(c*d+K+4*N*N+E)*R
    linear=psi*(E*T+pressure+curl)
    transport=2*(dr+d+r)
    wave=transport*U*U
    mean=(transport+2*N)*M*U
    return dict(linear=linear,complete_wave_wave=wave,complete_mean_wave=mean,
                F=linear+wave+mean,E=E,pressure_gradient=psi*pressure,
                principal_on_curl=psi*(c*d+K+4*N*N)*R,input_order=n+3)


class PrefixBudget:
    """Same real signed wave; fresh fast then fresh slow maps, never a reset."""
    def __init__(self,fixed,initial,wave,linear,KT1,KA1,time_inverse,slow,trace=None):
        self.fixed,self.initial,self.wave,self.linear=fixed,initial,wave,linear
        self.KT1,self.KA1,self.time_inverse,self.slow=KT1,KA1,time_inverse,slow
        self.trace=trace

    def data(self,name,n):
        order(n)
        if self.trace is not None:self.trace.add((name,n))
        return getattr(self,name)(n)

    @lru_cache(None)
    def rows(self,n,after_fast=False):
        state=self.postfast(n+3) if after_fast else self.data("initial",n+3)
        wave=self.data("wave",n+3)
        return full_mean(n,dict(self.data("fixed",n+3),**state,w=wave["U"]*wave["U"]))

    @lru_cache(None)
    def fast(self,n):
        order(n);time=self.data("time_inverse",n);a=self.data("KA1",n)
        coeff=comb(n+5,5)*24*2**(n+4)*time
        tang=coeff*factorial(n+4)*(self.rows(n+4)["Etheta"]+(1+a)*self.rows(n+4)["Ez"])
        rad=comb(n+5,5)*self.data("KT1",n+1)*24*2**(n+5)*time*factorial(n+5)*self.rows(n+5)["Ez"]
        return dict(x=rad,u=tang,z=tang)

    @lru_cache(None)
    def postfast(self,n):
        old=self.data("initial",n);inc=self.fast(n)
        return {k:old[k]+inc[k] for k in STATE}

    @lru_cache(None)
    def slow_increment(self,n):
        order(n);a=self.rows(n,True);anext=self.rows(n+1,True)
        k=self.data("slow",n);knext=self.data("slow",n+1)
        return dict(u=k["H"]*(k["Kp"]*a["P"]+k["Kz"]*a["Jz"]),
                    z=k["H"]*k["Kg"]*a["Jtheta"],
                    x=(n+1)*knext["H"]*knext["Kpsi"]*anext["Jtheta"])

    @lru_cache(None)
    def final(self,n):
        old=self.postfast(n);inc=self.slow_increment(n)
        return {k:old[k]+inc[k] for k in STATE}

    @lru_cache(None)
    def source(self,n):
        order(n);N=n+3
        fixed=self.data("fixed",N);wave=self.data("wave",N);state=self.final(N)
        data={k:fixed[k] for k in ("b","G","graph","c","r")}
        return source_recipe(n,dict(data,**self.data("linear",N),**wave,M=sum(state.values(),P())))

    def raw_source(self,n):
        order(n)
        return factorial(n)*self.source(n)["F"]

    def inverse_source(self,inverse_module):
        # Independently loaded historical helpers have distinct P classes.
        # Transfer exact coefficient dictionaries, with no numerical coercion.
        return lambda n:inverse_module.P(self.raw_source(n).c)


def fixture(trace=None):
    return PrefixBudget(fixed=lambda n:{k:P({0:i+1,1:Q(n+1,5)}) for i,k in enumerate(FIXED)},
        initial=lambda n:{k:P({0:i+1,1:Q(n+1,7)}) for i,k in enumerate(STATE)},
        wave=lambda n:{k:P({0:i+1,1:Q(n+1,11)}) for i,k in enumerate(("T","R","Pi","U"))},
        linear=lambda n:{k:P({0:i+1,1:Q(n+1,13)}) for i,k in enumerate(LINEAR)},
        KT1=lambda n:P(n+2),KA1=lambda n:P(n+1),time_inverse=lambda n:6*P({1:1}),
        slow=lambda n:{k:P({0:i+1,1:Q(n+1,17)}) for i,k in enumerate(("H","Kp","Kz","Kg","Kpsi"))},trace=trace)
