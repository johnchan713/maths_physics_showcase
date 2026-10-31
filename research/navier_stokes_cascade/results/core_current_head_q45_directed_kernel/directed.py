"""Fixed physical plateau witnesses -> same-pole Q44 integral constants."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
radial=load('q45_radial','core_current_head_q44_radial_multipliers/radial.py')
fixed=load('q45_fixed','core_current_head_q30_fixed_support_join/fixed_maps.py')
EP=radial.provider.EP
ExpSum=radial.provider.fast.refresh.fresh.ExpSum
def exponential(t):return EP(ExpSum([(1,t)]))

class Shell:
    def __init__(self,Xa,Xb,c0,c1,width):
        Xa,Xb,c0,c1,width=map(Q,(Xa,Xb,c0,c1,width))
        self.gap=fixed.gap(Xa,Xb,c0,c1)
        if not 0<width<=1:raise ValueError('Original positive flat width required')
        self.a,self.b=width**2,Q(4)
        self.lo=2*(Xb-Xa)/(Xb+Xa);self.hi=(Xb-Xa)/Xa
    def collar(self,p):
        radial.order(p)
        return min(Q(1),self.lo/4,self.gap,self.a/(1+p+self.a+self.b))
    @lru_cache(None)
    def bounds(self,p):
        r=self.collar(p);a,b=self.a,self.b
        E=exponential(Q(p*p,4)/min(a,b))
        K=r*(exponential(b/(self.lo-r)**2)+exponential(a/(self.lo-r)**2))
        K+=self.hi*E*exponential((a+b)/r**2)
        return dict(length=self.hi,edge=E,kernel=K,interior=exponential((a+b)/self.gap**2))
    def remaining(self,cutoff):
        def at(n,p):
            v=cutoff(n,p)
            if set(v)!={'order','pole','cutoff_order','bounds'} or set(v['bounds'])!={'chi','one_minus_chi','chi_U'}:
                raise ValueError('Only the original three raw cutoff norms are accepted')
            return dict(v,bounds=dict(v['bounds'],**self.bounds(p)))
        return at

def record():
    return dict(result='fixed-plateau-directed-kernel',uniform_gap_from_Q30=True,
                derived=['length','edge','kernel','interior'],cutoff_seminorms_required=True,
                A1_tail_actual_pulse_and_complete_residual_certified=False)
