"""Exact rational hull and raw radial multiplier recipes for Q42/Q43."""
from fractions import Fraction as Q
from math import ceil
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q44_chart',HERE.parent/'core_current_head_q43_compact_chart/chart.py')
chart=importlib.util.module_from_spec(s);s.loader.exec_module(chart)
provider=chart.provider
order=chart.order
REMAINING=provider.KEYS-{'f','g','k','jacobian'}

class RadialGeometry:
    def __init__(self,d,Rmin,Rmax):
        self.d,self.Rmin,self.Rmax=map(Q,(d,Rmin,Rmax))
        if not Q(1,3)<self.d<Q(5,2) or not 0<self.Rmin<=self.Rmax:raise ValueError('Original exponent and positive R hull required')
        self.Umin=min(Q(1),self.Rmin**3)
        self.Umax=max(Q(1),self.Rmax**3)
    def amplitude(self,b):
        m=ceil(abs(Q(b)))
        return max(Q(1),self.Rmin**(-m),self.Rmax**m)
    def power(self,n,b,scale=1):
        """Raw through-n U derivatives of scale*R**b, U=R**d."""
        order(n);b,scale=Q(b),Q(scale)
        if scale<0:raise ValueError('Nonnegative multiplier scale required')
        a=b/self.d;fall=Q(1);v=[Q(1)]
        for j in range(1,n+1):
            fall*=a-(j-1);v.append(abs(fall)*self.Umin**(-j))
        return scale*self.amplitude(b)*max(v)
    def multipliers(self,n,e=1):
        order(n)
        if type(e) is not int or e not in (0,1,2):raise ValueError('Original radial moment channel required')
        d=self.d
        return dict(f=self.power(n,e-d+1,1/d),g=self.power(n,-e),
                    k=self.power(n,d-1-e,d),jacobian=d*self.Umax/2)

def build(d,Rmin,Rmax,remaining,time_inverse,max_order,max_pole,*,e=1):
    """Same hull feeds both raw multipliers and the chart roundtrip.

    remaining(n,p) supplies exactly Q42's declared coverage and seven
    non-radial bounds. It retains actual cutoff/kernel/weight premises.
    e=1 is the axial fast map; e=0/2 are separate moment channels.
    """
    g=RadialGeometry(d,Rmin,Rmax);g.multipliers(0,e)
    def raw(n,p):
        a=remaining(n,p)
        if set(a)!={'order','pole','cutoff_order','bounds'} or set(a['bounds'])!=REMAINING:raise ValueError('Only original cutoff/kernel/weight inputs are accepted')
        return dict(a,bounds=dict(a['bounds'],**g.multipliers(n,e)))
    c=provider.CompactProvider(raw,time_inverse,max_order,max_pole)
    return chart.ChartedCompact(c,d=g.d,Rmin=g.Rmin,Rmax=g.Rmax,Umin=g.Umin,Umax=g.Umax)

def record():
    return dict(result='derived-radial-hull-and-compact-multipliers',
        Umin='min(1,Rmin^3)',Umax='max(1,Rmax^3)',derived_inputs=['f','g','k','jacobian'],
        chart_and_multiplier_hull_identical=True,cutoff_kernel_and_time_inputs_retained=True,
        initial_A1_tail_full_source_residual_or_blowup_certified=False)
