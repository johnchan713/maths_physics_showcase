"""Native raw cutoff jets from one physical X profile, with slow scale motion."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
from functools import lru_cache
HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
directed=load('q46_directed','core_current_head_q45_directed_kernel/directed.py')
coordinates=load('q46_coordinates','core_primary_target_cutoff_budget/bounds.py')
radial=directed.radial

class Moving:
    def __init__(self,d,Rmin,Rmax,Xb,width,shape):
        self.geo=radial.RadialGeometry(d,Rmin,Rmax)
        self.coords=coordinates.ChartBudget(Rmin,width)
        self.Xb=Q(Xb);self.shape=tuple(map(Q,shape))
        if self.Xb<=0 or not self.shape or min(self.shape)<0:raise ValueError('Positive shell end and nonnegative normalized shape seminorms required')
    @lru_cache(None)
    def raw(self,n):
        radial.order(n)
        if n>=len(self.shape):raise ValueError('Original shape derivative ledger too short')
        # Exact-degree sums of normalized jets, then one raw conversion.
        L=[0]+[2/(self.geo.d*j*self.geo.Umin**j)+self.coords.log_s(j) for j in range(1,n+1)]
        X=[self.Xb]
        for k in range(1,n+1):
            X.append(sum(Q(j,k)*L[j]*X[k-j] for j in range(1,k+1)))
        power=[1]+[0]*n;out=[self.shape[0]]+[0]*n
        for k in range(1,n+1):
            power=[sum(X[j]*power[i-j] for j in range(1,i+1)) for i in range(n+1)]
            for i in range(k,n+1):out[i]+=self.shape[k]*power[i]
        return tuple(factorial(i)*v for i,v in enumerate(out))
    def __call__(self,n,p):
        radial.order(n);radial.order(p);raw=self.raw(n+1)
        return dict(order=n,pole=p,cutoff_order=n+1,bounds=dict(
            chi=max(raw[:n+1]),one_minus_chi=max(1+raw[0],*raw[1:n+1]) if n else 1+raw[0],
            chi_U=max(raw[1:])))
def build(*,d,Rmin,Rmax,Xa,Xb,c0,c1,width,shape,time_inverse,max_order,max_pole,e=1):
    shell=directed.Shell(Xa,Xb,c0,c1,width)
    jets=Moving(d,Rmin,Rmax,Xb,width,shape)
    if len(jets.shape)<max_order+2:raise ValueError('Reserve the additional cutoff derivative')
    return radial.build(d,Rmin,Rmax,shell.remaining(jets),time_inverse,max_order,max_pole,e=e)

def record():
    return dict(result='moving-physical-cutoff-jets',fixed_shape_ceiling='compact order + 1',
        shifted_source_undifferentiated_in_M=True,scale_derivatives_retained=True,
        actual_shape_values_initial_A1_tail_or_full_residual_certified=False)
