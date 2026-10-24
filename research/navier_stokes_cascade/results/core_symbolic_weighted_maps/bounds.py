"""Finite positive exponential constants for the unchanged weighted maps."""
from fractions import Fraction as Q
from math import ceil,comb,factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('symbolic_map_rows',HERE.parent/'core_initial_haar_row_budget/budgets.py')
rows=importlib.util.module_from_spec(sp);sp.loader.exec_module(rows)
exact,order=rows.exact,rows.order

class ExpSum:
    """Exact sum c*exp(e), rational c,e>=0. No numerical exponential."""
    def __init__(self,terms):
        self.terms={}
        for c,e in terms:
            c,e=exact(c),exact(e)
            if min(c,e)<0:raise ValueError('Nonnegative sufficient bound required')
            if c:self.terms[e]=self.terms.get(e,Q(0))+c
    def __add__(self,other):return ExpSum([(c,e) for obj in (self,other) for e,c in obj.terms.items()])
    def __mul__(self,other):return ExpSum([(c*d,e+f) for e,c in self.terms.items() for f,d in other.terms.items()])
    def scale(self,c):return self*ExpSum([(c,0)])
    def rational_upper(self,max_power=10000):
        order(max_power)
        if any(ceil(e)>max_power for e in self.terms):raise ValueError('Keep giant fixed constants symbolic')
        return sum((c*3**ceil(e) for e,c in self.terms.items()),Q(0))
    def record(self):return [[str(c),str(e)] for e,c in sorted(self.terms.items())]

def edge_absorption(p,width):
    p,width=map(exact,(p,width))
    if p<0 or not 0<width<=1:raise ValueError('Original width and nonnegative pole required')
    N=ceil(p/2)
    return max(Q(1),Q(factorial(N))*(2/width**2)**N)

def weighted_maps(p,L,width,gap,*,jacobian,KT_profile,KA_profile):
    """Profiles are the original raw multiplier/cutoff product constants.

    KT_profile=2^(3m) N_m(g)[N_m(1-chi)+N_m(chi)]N_m(f).
    KA_profile is the original UNWEIGHTED A_e constant (including L).
    jacobian bounds dU/dx on the same whole shell. All inputs fixed.
    """
    p,L,width,gap,J,kt,ka=map(exact,(p,L,width,gap,jacobian,KT_profile,KA_profile))
    if p<0 or L<=0 or not 0<width<=1 or not 0<gap<=L/2 or min(J,kt,ka)<0:raise ValueError('Positive shell, actual plateau gap and nonnegative fixed profile bounds required')
    r=min(Q(1),L/4,gap);a,b=width**2,Q(4);q=max(0,ceil((p-3)/2))
    def D(c):return sum((Q(comb(q,j)*factorial(j))*(r*r/c)**j for j in range(q+1)),Q(0))/(2*c)
    E=edge_absorption(p,width)
    # Sum safely dominates the historical max of three branches.
    C=ExpSum([(D(a),b/(L-r)**2),(D(b),a/(L-r)**2),(L*E,(a+b)/(r*r))])
    Z=ExpSum([(1,(a+b)/(gap*gap))])
    return dict(KT=C.scale(J*kt),KA=Z.scale(E*ka),primitive=C,
                proof_collar=r,edge_absorption=E,weight_floor=Z)

def fast_coefficients(m,KT,KA,*,Et_raw,Ez_raw,Ez_radial_raw,time_floor):
    """Raw source bounds at m+4/m+5; return normalized output coefficients.

    time_floor bounds c^-1/S. The S factor remains EXTERNAL and exact.
    KT is at m+1, KA at m, as in the historical fast velocity identity.
    """
    order(m)
    if not isinstance(KT,ExpSum) or not isinstance(KA,ExpSum):raise TypeError('Symbolic weighted operators required')
    et,ez,ezr,c=map(exact,(Et_raw,Ez_raw,Ez_radial_raw,time_floor))
    if min(et,ez,ezr,c)<0:raise ValueError('Nonnegative source and time-floor coefficients required')
    gain=comb(m+5,5)*24*2**(m+4)*c
    tangential=(ExpSum([(et+ez,0)])+KA.scale(ez)).scale(gain)
    radial=KT.scale(ezr*comb(m+5,5)*24*2**(m+5)*c)
    return dict(Ft=tangential,Fr=radial,external_S_power=1,
                raw_source_orders=(m+4,m+5),actual_axial_A1_retained=True)
