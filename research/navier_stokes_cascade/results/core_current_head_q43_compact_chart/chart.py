"""Finite same-weight radial chart conversion for Q42 compact operators."""
from fractions import Fraction as Q
from math import factorial
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q43_provider',HERE.parent/'core_current_head_q42_compact_provider/provider.py')
provider=importlib.util.module_from_spec(s);s.loader.exec_module(provider)
order=provider.order

def binomial(a,j):
    b=Q(1)
    for k in range(j):b*=Q(a-k,k+1)
    return b

def transfer(n,exponent,source_max,target_min,normalized_input):
    """source coordinate = target coordinate ** exponent, four identities.

    source_max and target_min enclose the actual common positive hull.
    Return a sufficient raw-output factor (also bounds normalized output).
    """
    order(n);a,H,l=map(Q,(exponent,source_max,target_min))
    if min(a,H,l)<=0:raise ValueError('Positive chart and hull bounds required')
    if type(normalized_input) is not bool:raise TypeError('Declare the input norm')
    jets=[Q(0)]+[H*abs(binomial(a,j))/l**j+4*int(j==1) for j in range(1,n+1)]
    e=[Q(1)]
    for k in range(1,n+1):
        e.append(sum((jets[j]*e[k-j]*(1 if normalized_input else Q(j,k)) for j in range(1,k+1)),Q(0)))
    return factorial(n)*sum(e)

class ChartedCompact:
    """Normalized R -> raw U -> Q42 operator -> normalized R.

    Enclosures must describe the same fixed hull U=R**d. No new geometry
    or source support is chosen. The time inverse input remains in the
    original Q40 convention; only KT/KA are converted here.
    """
    def __init__(self,raw,*,d,Rmin,Rmax,Umin,Umax):
        if not isinstance(raw,provider.CompactProvider):raise TypeError('Q42 raw compact provider required')
        self.raw=raw;self.d,self.Rmin,self.Rmax,self.Umin,self.Umax=map(Q,(d,Rmin,Rmax,Umin,Umax))
        if not Q(1,3)<self.d<Q(5,2) or not 0<self.Rmin<=self.Rmax or not 0<self.Umin<=self.Umax:raise ValueError('Original exponent and positive hull required')
    @lru_cache(None)
    def factors(self,n):
        return (transfer(n,1/self.d,self.Rmax,self.Umin,True),
                transfer(n,self.d,self.Umax,self.Rmin,False))
    @lru_cache(None)
    def __call__(self,n,pole):
        v=self.raw(n,pole);a,b=self.factors(n)
        return dict(v,KT=a*b*v['KT'],KA=a*b*v['KA'])

def record():
    return dict(result='same-weight-compact-chart-roundtrip',forward='n! sum [t^j](1-A(t))^-1',
        reverse='n! sum [t^j]exp(B(t))',same_pole=True,extra_epsilon_loss=0,
        derivative_order_loss=0,actual_enclosures_profile_jets_and_A1_tail_required=True,
        complete_source_residual_or_blowup_certified=False)
