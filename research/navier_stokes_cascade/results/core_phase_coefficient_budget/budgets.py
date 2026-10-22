"""Actual frozen-phase coefficient recipes from full background jet bounds."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module);return module


inverse=load("phase_budget_inverse",HERE.parent/"core_first_sourced_inverse_budget/budgets.py")
P,order,exact=inverse.P,inverse.order,inverse.exact


def product(a,b,n):
    order(n)
    return sum((comb(n,j)*a(j)*b(n-j) for j in range(n+1)),P())


def carrier_bounds(lambda_max,u,gmin,Rmax,Lmin):
    lam,u,g,R,L=map(exact,(lambda_max,u,gmin,Rmax,Lmin))
    if min(lam,g,R,L)<=0 or u<0:raise ValueError("Positive growth/shear/radius/length and nonnegative tilt bounds required")
    B=max(Q(1),lam)  # B >= sqrt(lambda_max); no square root computation
    qz=B*(1+u/(L*g))
    return dict(theta=R*qz+1,axial=qz,radial=B*u/2)


def normal_gate(lambda_min,u,Cphi,S,epsilon_sum):
    lam,u,C,S,E=map(exact,(lambda_min,u,Cphi,S,epsilon_sum))
    if lam<=0 or min(u,C,E)<0 or S<1:raise ValueError("Positive growth floor and valid value-comparison inputs required")
    beta=lam/(4*(1+u*u)**2)
    gates=dict(whole_path_rounding_and_axial_error=S*S*E<=1,
               full_normal_comparison_margin=36*C*C<=beta*S*S)
    return dict(gates=gates,beta=beta,inverse_normal_square=4/beta,phase_error=3*C/S,
                tilt_upper=u,lambda_min=lam,
                assumptions="actual comparison (5), epsilon_sum=epsilon+epsilon^2+k^-1; no new accepted band selected")


class PhaseBudget:
    def __init__(self,F,G,btilde,path,Rmin,Rmax,Lmin,u,lambda_max,carriers,trace=None):
        self.F,self.G,self.btilde,self.path=F,G,btilde,path
        self.Rmin,self.Rmax,self.Lmin,self.u,self.lam=map(exact,(Rmin,Rmax,Lmin,u,lambda_max))
        if min(self.Rmin,self.Rmax,self.Lmin,self.lam)<=0 or self.Rmin>self.Rmax or self.u<0:
            raise ValueError("Valid positive shell, path and growth bounds required")
        if set(carriers)!={"theta","axial","radial"}:raise ValueError("All three frozen carrier bounds required")
        self.p,self.pz,self.x0=(exact(carriers[k]) for k in ("theta","axial","radial"))
        if min(self.p,self.pz,self.x0)<0:raise ValueError("Nonnegative absolute carrier bounds required")
        self.trace=trace

    def data(self,name,n):
        order(n)
        if self.trace is not None:self.trace.add((name,n))
        v=getattr(self,name)(n)
        if not isinstance(v,P):raise TypeError("Positive full-background raw mixed-jet polynomials required")
        return v

    def radius(self,n):
        order(n);return P(self.Rmax if n==0 else int(n==1))

    def reciprocal_radius(self,n):
        order(n);return P(Q(factorial(n),1)/self.Rmin**(n+1))

    @lru_cache(None)
    def Hder(self,n):
        order(n)
        return self.p*self.data("F",n+1)+self.pz*self.data("G",n+1)

    @lru_cache(None)
    def components(self,n):
        order(n);vh=product(lambda j:self.data("path",j),self.Hder,n)
        return dict(radial=int(n==0)*P(self.x0)+vh,
                    angular=self.p*self.reciprocal_radius(n),
                    axial=int(n==0)*P(self.pz)+vh)  # epsilon<=1

    def normal(self,n):return sum(self.components(n).values(),P())
    def moving_normal(self,n):return 2*self.Hder(n)

    @lru_cache(None)
    def matrix(self,n):
        order(n)
        return 4*self.data("F",n)+product(self.radius,lambda j:self.data("F",j+1),n)+self.data("G",n+1)

    @lru_cache(None)
    def eikonal(self,n):
        order(n)
        transport=lambda j:self.Hder(j)+product(lambda k:self.data("G",k),self.Hder,j)
        return product(lambda j:self.data("path",j),transport,n)+product(
            lambda j:self.data("btilde",j),lambda j:self.components(j)["radial"],n)

    def bR(self,n):return self.data("btilde",n+1)

    def Bz(self,n):
        order(n)
        return self.data("btilde",n+1)+product(self.radius,lambda j:self.data("F",j+1),n)+self.data("G",n+1)

    def sjet(self,n):
        order(n)
        return P(Q(3,2)*self.u) if not n else self.u/self.Lmin*self.data("path",n)

    @lru_cache(None)
    def shape(self,n):
        return P(int(n==0))+product(self.sjet,self.sjet,n)

    @lru_cache(None)
    def increment_power(self,n,k):
        if not k:return P(int(n==0))
        if n<k:return P()
        return sum((Q(1,factorial(j))*self.shape(j)*self.increment_power(n-j,k-1)
                    for j in range(1,n+1)),P())

    @lru_cache(None)
    def growth_shape(self,n):
        order(n)
        if not n:return P(1)
        return factorial(n)*sum((Q(comb(2*k,k),4**k)*self.increment_power(n,k)
                                 for k in range(1,n+1)),P())

    def profile_rate(self,n):
        # lambda(v)=lambda0/sqrt(1+s_sigma(v)^2), NOT lambda0.
        return self.lam*(self.growth_shape(n)+self.shape(n))

    def normalized(self,kind,n,dimension=5):
        order(n);order(dimension)
        if not dimension:raise ValueError("At least one coefficient variable required")
        callback=getattr(self,kind)
        return sum((Q(dimension**j,factorial(j))*callback(j) for j in range(n+1)),P())

    def source_coefficients(self,n):
        return dict(b=self.normalized("b_raw",n),G=self.normalized("G_raw",n),
                    V=self.normalized("V_raw",n),r=self.normalized("reciprocal_radius",n),
                    normal=self.normalized("normal",n),K=self.normalized("matrix",n),
                    Eik=self.normalized("eikonal",n),bR=self.normalized("bR",n),Bz=self.normalized("Bz",n))

    def b_raw(self,n):return self.data("btilde",n)
    def G_raw(self,n):return self.data("G",n)
    def V_raw(self,n):return product(self.radius,lambda j:self.data("F",j),n)

    def forced(self,source,tangent,length,gate,harmonic=1,module=inverse):
        if not all(gate["gates"].values()):raise ValueError("Actual full-normal value gate is not established")
        if gate["tilt_upper"]<self.u or gate["lambda_min"]>self.lam:
            raise ValueError("Value gate and coefficient growth/tilt inputs are inconsistent")
        cv=lambda fn:lambda n:module.P(fn(n).c)
        return module.ForcedBudget(normal=cv(self.normal),moving_normal=cv(self.moving_normal),
            matrix=cv(self.matrix),source=cv(source),profile_rate=cv(self.profile_rate),
            tangent=module.P(tangent.c),length=module.P(length.c),
            inverse_normal_square=gate["inverse_normal_square"],harmonic=harmonic)


def fixture(trace=None):
    return PhaseBudget(F=lambda n:P(n+2),G=lambda n:P(n+3),btilde=lambda n:P(n+1),
        path=lambda n:P({1:2}) if n==0 else P(3 if n==1 else 0),
        Rmin=1,Rmax=2,Lmin=1,u=2,lambda_max=3,
        carriers=carrier_bounds(3,2,Q(1,2),2,1),trace=trace)
