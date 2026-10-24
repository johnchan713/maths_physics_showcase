"""Same-weight raw U jets to factorial R jets on the original radial hull."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('reverse_chart_bridge',HERE.parent/'core_weighted_primary_covariance_join/bounds.py')
bridge=importlib.util.module_from_spec(sp);sp.loader.exec_module(bridge)
P,order,exact=bridge.P,bridge.order,bridge.primary.exact

class ReverseChart:
    """U=R**d_r, 1/3<d_r<5/2, with four other coordinates unchanged.

    Input: maximum RAW derivative through N at the corresponding U point.
    Output: both the sum of normalized derivatives and the sum of RAW
    derivatives through N in R. Same physical weight and pole, no quotient
    differentiation. All bounds are finite constants independent of S/epsilon.
    """
    def __init__(self,Rmin,Rmax):
        self.Rmin,self.Rmax=map(exact,(Rmin,Rmax))
        if not 0<self.Rmin<=self.Rmax:raise ValueError('Same closed positive radial hull required')
        self.B=max(Q(1),self.Rmax**2,1/self.Rmin)

    def radial(self,j):
        """Normalized positive-order jet of U(R), from the original g bound."""
        order(j)
        if not j:raise ValueError('Only positive-order coordinate increments are used')
        return 3*self.B/self.Rmin**(j-1)

    def increment(self,j):
        order(j)
        return Q(0) if not j else 4*int(j==1)+self.radial(j)

    @lru_cache(None)
    def exponential(self,n):
        """[t**n] exp(a(t)), a=4t+sum radial(j)t**j; exact rational."""
        order(n)
        if not n:return Q(1)
        return sum((Q(j,n)*self.increment(j)*self.exponential(n-j) for j in range(1,n+1)),Q(0))

    def transfer(self,N):
        order(N)
        return P(factorial(N)*sum((self.exponential(n) for n in range(N+1)),Q(0)))

class ChartedPrimary:
    """Discharge the bridge's tensor chart slot, retain every other R input."""
    def __init__(self,primary,chart):
        if not isinstance(primary,bridge.WeightedPrimary) or not isinstance(chart,ReverseChart):raise TypeError('Same-family primary recipe and explicit original-hull chart required')
        self.primary,self.chart=primary,chart

    def covariance_R(self,n):
        c=self.primary.covariance(n);t=self.chart.transfer(n)
        return dict(c,W=t*c['W'],Dcurl=t*c['Dcurl'])

    def initial_rows(self,n,data_R,higher_stress_R):
        """Other fields/operators must already obey the historical R contracts."""
        order(n)
        result=self.primary.initial_rows(n,data_R,higher_stress_R,self.chart.transfer(n+3))
        result.update(tensor_chart_transfer_derived=True,other_R_inputs_still_required=True)
        return result

def from_original(p,selector_inputs,ell,**kw):
    """Take the hull from the SAME returned complete-curl operator, not a new hull."""
    bundle,primary=bridge.from_original(p,selector_inputs,ell,**kw)
    operator=bundle[2]
    return bundle,ChartedPrimary(primary,ReverseChart(operator.Rmin,operator.Rmax))
