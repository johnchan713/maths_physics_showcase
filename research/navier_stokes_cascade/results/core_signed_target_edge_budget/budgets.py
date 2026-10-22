"""Fixed-width two-edge constants and current projected-target norm routing."""
from fractions import Fraction as Q
from functools import lru_cache
from math import ceil,comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("signed_edge_root",HERE.parent/"core_primary_root_budget/budgets.py")
root=importlib.util.module_from_spec(spec);spec.loader.exec_module(root)
P=root.P;atlas=root.atlas;exact=atlas.exact;order=atlas.order


def primitive_recipe(p,L,width):
    """Return three coeff*exp(exponent) branches; never evaluate giant widths."""
    p,L,width=map(exact,(p,L,width))
    if p<0 or L<=0 or not 0<width<=1:raise ValueError("Nonnegative pole, positive shell length and width in (0,1] required")
    r=min(Q(1),L/4);a,b=width**2,Q(4);q=max(0,ceil((p-3)/2))
    def D(c):return sum((Q(comb(q,j)*factorial(j))*(r*r/c)**j for j in range(q+1)),Q(0))/(2*c)
    return dict(pole=str(p),L=str(L),width=str(width),collar=str(r),kernel_degree=q,
                branches=dict(left=(D(a),b/(L-r)**2),right=(D(b),a/(L-r)**2),
                              middle=(L*atlas.edge_constant(p,width),(a+b)/(r*r))))


def rational_upper(recipe,max_power=10000):
    """Diagnostic materialization only: e^x<=3^ceil(x), x>=0."""
    order(max_power)
    branches=recipe["branches"]
    if any(ceil(exponent)>max_power for coefficient,exponent in branches.values()):
        raise ValueError("Keep the fixed giant-width constant symbolic")
    return max(coefficient*3**ceil(exponent) for coefficient,exponent in branches.values())


class TargetBudget:
    """All inputs are sufficient positive-polynomial bounds, not actual fields."""
    def __init__(self,row,bump,moment,primitive,radial,reciprocal,poles,chart=None,trace=None):
        self.row,self.bump,self.moment=row,bump,moment
        self.primitive,self.radial,self.reciprocal=primitive,radial,reciprocal
        self.poles=poles;self.trace=trace
        self.chart=chart or (lambda n:P(int(n==1)))

    def data(self,name,n):
        order(n)
        if self.trace is not None:self.trace.add((name,n))
        value=getattr(self,name)(n)
        if not isinstance(value,P):raise TypeError("Positive scale-polynomial bounds required")
        return value

    @lru_cache(None)
    def H(self,n):
        order(n)
        return self.data("row",n)+sum((self.data("bump",j)*self.data("moment",n-j) for j in range(n+1)),P())

    @lru_cache(None)
    def F(self,n):
        order(n)
        return self.data("primitive",n)*self.H(n)+(self.H(n-1) if n else P())

    @lru_cache(None)
    def K(self,n):
        order(n)
        return sum((self.data("radial",j)*self.F(n-j) for j in range(n+1)),P())

    @lru_cache(None)
    def Vx(self,n):
        order(n)
        return sum((self.K(j)*self.data("reciprocal",n-j) for j in range(n+1)),P())

    @lru_cache(None)
    def chart_power(self,n,k):
        order(n);order(k)
        if k==0:return P(int(n==0))
        if n<k:return P()
        return sum((self.data("chart",j)*self.chart_power(n-j,k-1) for j in range(1,n+1)),P())

    @lru_cache(None)
    def V(self,n):
        order(n)
        return sum((self.Vx(k)*self.chart_power(n,k) for k in range(n+1)),P())

    def edge(self,N):
        order(N);p=[exact(self.poles(n)) for n in range(N+1)]
        if any(x<0 for x in p) or any(p[n]<p[n-1] for n in range(1,N+1)):
            raise ValueError("Nonnegative nondecreasing incoming edge poles required")
        r=[max(p[j]+3*(n-j) for j in range(n+1)) for n in range(N+1)]
        return dict(source=p,quotient=r,d=max([Q(0)]+[r[n]-3*n for n in range(N+1)]))


def fixture(trace=None):
    return TargetBudget(row=lambda n:P({0:n+1,1:2}),bump=lambda n:P(n+2),
                        moment=lambda n:P({1:n+3}),primitive=lambda n:P(rational_upper(primitive_recipe(n+3,4,Q(1,2)))),
                        radial=lambda n:P(n+1),reciprocal=lambda n:P(2**n),poles=lambda n:6*n+3,
                        chart=lambda n:P(n+1),trace=trace)
