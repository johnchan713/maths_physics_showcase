"""Positive source coefficient compression; native leaf estimates remain premises."""
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q113_compact_shape_scope/q113.py'
s=importlib.util.spec_from_file_location('q114_scope',p);scope=importlib.util.module_from_spec(s);s.loader.exec_module(scope)
radial=scope.profiles.charts.radial
order=scope.profiles.order

@dataclass(frozen=True)
class Growth:
    """P(S)<=exp(height)*(1+S)^degree for S>=0; None degree is exact zero."""
    degree: object
    height: Q=Q(0)
    def __post_init__(self):
        if self.degree is not None:order(self.degree)
        object.__setattr__(self,'height',Q(self.height))
        if self.height<0:raise ValueError('Nonnegative fixed log-height required')
    def __add__(self,other):
        other=as_growth(other)
        if self.degree is None:return other
        if other.degree is None:return self
        return Growth(max(self.degree,other.degree),1+max(self.height,other.height))
    __radd__=__add__
    def __mul__(self,other):
        other=as_growth(other)
        if self.degree is None or other.degree is None:return Growth(None)
        return Growth(self.degree+other.degree,self.height+other.height)
    __rmul__=__mul__

def as_growth(value):
    if isinstance(value,Growth):return value
    value=Q(value)
    if value<0:raise ValueError('Positive coefficient algebra only')
    return Growth(0,value) if value else Growth(None)

def compress(poly):
    """Compress an existing EP, preserving exp(e) as a FIXED coefficient."""
    poly=radial.EP(poly);out=Growth(None)
    for degree,coeff in poly.c.items():
        for exponent,c in coeff.terms.items():
            out+=Growth(degree,Q(c)+Q(exponent))
    return out

def reduced(n,mean,W,a,fixed,absorption):
    """Same positive reduced polynomial as Q104, in the Growth semiring."""
    order(n);a=Q(a)
    if a<=0 or set(mean)!={'x','u','z'} or set(fixed)!={'b','V','G','c','graph'}:
        raise ValueError('Original positive hull and complete source leaves required')
    x,u,z=(as_growth(absorption)*as_growth(mean[k]) for k in ('x','u','z'))
    b,V,G,c=(as_growth(fixed[k]) for k in ('b','V','G','c'))
    W=as_growth(W);d=n+1;dd=(n+1)*(n+2)
    r=sum((a**(-j-1) for j in range(n+1)),Q(0))
    S=2*b*x+x*x+W;H=b*z+G*x+x*z+W;C=2*V*u+u*u+W
    return d*(1+4*c)*x+r*S+d*H+r*C+dd*x

class SourceGrowth:
    def __init__(self,kernels):
        if not isinstance(kernels,scope.CompactShapeScope):raise TypeError('Q113 same-state kernel scope required')
        self.kernels=kernels;self.state=kernels.state
    def source(self,J,bundle):
        self.kernels.plan(0) # recheck the unchanged current state
        out=self.kernels.kernels.bridge.radial.source(J,bundle)
        if out['source_state'] is not self.state:raise ValueError('Current source state changed')
        return dict(growth=compress(out['coefficient']),source=out,
                    native_leaf_inequalities_still_required=True,actual_A0_gain_certified=False)

def record():
    return dict(status='current-reduced-source-polynomial-degree-and-height',
                bound='exp(H_J)*(1+S)^D_J',
                fixed_exponential_coefficients_not_exp_in_S=True,
                reduced_source_graph_degree_absent=True,
                source_native_contracts_or_A0_gain_certified=False,
                full_cone_pulse_residual_or_blowup=False)
