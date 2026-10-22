"""Positive fixed-order root/edge/curl budgets; inputs are bounds, not fields."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("root_budget_atlas",HERE.parent/"core_primary_aggregate_atlas/proof_checks.py")
atlas=importlib.util.module_from_spec(spec);spec.loader.exec_module(atlas)
P=atlas.P;order=atlas.order;exact=atlas.exact


def convolution(a,b,n):
    return sum((a(j)*b(n-j) for j in range(n+1)),P())


class RootBudget:
    """All jets are total-degree Taylor sums; raw norms are converted last."""
    def __init__(self,H,target,wave,curl,loghalf,K0,beta,R,mass,width,trace=None):
        self.H,self.target,self.wave,self.curl,self.loghalf=H,target,wave,curl,loghalf
        self.K0,self.beta,self.R,self.mass,self.width=map(exact,(K0,beta,R,mass,width))
        if min(self.K0,self.beta,self.R,self.mass)<=0 or self.R*self.R*self.beta<1:
            raise ValueError("Positive matrix/root/mass constants with R>=beta^-1/2 required")
        atlas.edge_constant(0,self.width)
        self.trace=trace

    def data(self,name,n):
        order(n)
        if self.trace is not None:self.trace.add((name,n))
        p=getattr(self,name)(n)
        if not isinstance(p,P):raise TypeError("Positive scale-polynomial input required")
        return p

    @lru_cache(None)
    def K(self,n):
        order(n)
        if not n:return P(self.K0)
        return self.K0*sum((self.data("H",j)*self.K(n-j) for j in range(1,n+1)),P())

    @lru_cache(None)
    def B(self,n):
        order(n)
        return sum((self.K(j)*self.data("target",n-j) for j in range(n+1)),P())

    @lru_cache(None)
    def power(self,n,k):
        if k==0:return P(int(n==0))
        if n<k:return P()
        return sum((self.B(j)*self.power(n-j,k-1) for j in range(1,n+1)),P())

    @lru_cache(None)
    def reciprocal(self,n):
        order(n)
        if not n:return P(self.R)
        coefficient=Q(1);result=P()
        for k in range(1,n+1):
            coefficient*=Q(2*k-1,2*k)
            result+=self.R*coefficient/self.beta**k*self.power(n,k)
        return result

    @lru_cache(None)
    def root(self,n):
        order(n)
        return convolution(self.B,self.reciprocal,n)

    @lru_cache(None)
    def E(self,n):
        order(n)
        if not n:return P(1)
        return sum((Q(j,n)*self.data("loghalf",j)*self.E(n-j) for j in range(1,n+1)),P())

    @lru_cache(None)
    def amplitude(self,n):
        order(n)
        return self.mass*P({1:1})*sum((self.E(i)*self.root(j)*self.data("wave",n-i-j)
                                      for i in range(n+1) for j in range(n-i+1)),P())

    @lru_cache(None)
    def complete_curl(self,n):
        order(n)
        return sum((self.data("curl",j)*(self.amplitude(n-j)+(n-j+1)*self.amplitude(n-j+1))
                    for j in range(n+1)),P())

    def native(self,n):
        order(n)
        a=2*factorial(n)*atlas.edge_constant(6*n+3,self.width)*sum((self.amplitude(j) for j in range(n+1)),P())
        b=2*factorial(n)*atlas.edge_constant(6*n+9,self.width)*sum((self.complete_curl(j) for j in range(n+1)),P())
        return a,b

    @lru_cache(None)
    def aggregate(self,n):
        a,b=self.native(n)
        # Already unweighted in native(): do not apply the edge constants twice.
        return atlas.aggregate(n,a,b,0,0,self.width)

    def signed(self,n,xi):
        order(n)
        def V(j):return sum((self.K(k)*xi(j-k) for k in range(j+1)),P())
        return self.mass*Q(1,2)*P({1:1})*sum((self.E(i)*V(j)*self.reciprocal(k)*self.data("wave",n-i-j-k)
                                             for i in range(n+1) for j in range(n-i+1)
                                             for k in range(n-i-j+1)),P())


def fixture(trace=None):
    return RootBudget(H=lambda n:P({0:n+1,1:1}),target=lambda n:P({0:n+2,1:2}),
                      wave=lambda n:P({0:1,1:n+1}),curl=lambda n:P({0:2,1:n+1}),
                      loghalf=lambda n:P(n+1),K0=2,beta=Q(1,4),R=2,mass=12,
                      width=Q(1,8),trace=trace)
