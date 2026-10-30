"""Same-profile radial maps at the actual Q37 pressure pole."""
from fractions import Fraction as Q
from math import ceil,comb,factorial
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('q38_update',HERE.parent/'core_current_head_q37_signed_covariance_pressure/signed_update.py')
update=importlib.util.module_from_spec(sp);sp.loader.exec_module(update)
EP,order=update.EP,update.order

class NewPoleMaps:
    """Bounded maps on smooth flat sources with exact projected moments.

    This is a norm composition, not verification of an actual source field.
    """
    def __init__(self,state):
        if not isinstance(state,update.SignedCovariance):raise TypeError('Q37 same-choice signed state required')
        self.state=state;self.target=state.target;self.profile=state.target.source.profile
        if self.profile.max_order<state.max_output+1:raise ValueError('Reserve the next pressure-map profile jet')
    @lru_cache(None)
    def kernel(self,p):
        order(p)
        lo,hi=self.target.Llo,self.target.Lhi
        r=min(Q(1),lo/4);a=self.target.width**2;b=Q(4)
        k=max(0,ceil(Q(p-3,2)))
        def D(c):return sum((Q(comb(k,j)*factorial(j))*(r*r/c)**j for j in range(k+1)),Q(0))/(2*c)
        return EP(update.signed.maps.ExpSum([(D(a),b/(lo-r)**2),(D(b),a/(lo-r)**2),
            (hi*update.signed.maps.edge_absorption(p,self.target.width),(a+b)/(r*r))]))
    def recipe(self,n,p):
        order(n);order(p);self.profile.check(n)
        a,b=self.profile.Rmin,self.profile.Rmax
        kp=1+b*Q(1,2)*self.kernel(p)
        def te(e):
            positive=sum(Q(comb(e,j))*b**(e-j) for j in range(min(e,n)+1))
            negative=sum(Q(comb(e+j-1,j))*a**(-e-j) for j in range(n+1))
            return positive*negative*kp
        bump=self.profile.floor*(self.profile.sigma(n,1)+self.profile.sigma(n,2))
        return dict(order=n,pole=p,Kp=kp,primitive=te(1)+te(2),bump=bump)
    def pressure(self,n):return self.state.pressure(n,self.recipe)

def record():
    return dict(result='same-profile-new-pole-radial-map-composition',
        new_pressure_pole_used=True,independent_Kp_primitive_bump_callbacks_removed=True,
        exact_zero_moments_and_flat_extensions_required=True,
        kernel_applied_before_differentiation=False,radial_jacobian_retained=True,
        profile_floor_once_per_bump=True,whole_pressure_memory_retained=True,
        actual_source_and_fresh_means_certified=False,complete_residual_or_blowup=False)
