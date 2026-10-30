"""Fixed-profile Q35 source to the original signed root; conditional bounds."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
pipeline=load('q36_pipeline','core_current_head_q35_initial_projected_source/pipeline.py')
signed=load('q36_signed','core_symbolic_signed_target_join/bounds.py')
EP,order,exact=signed.EP,signed.order,signed.exact

def rebind(value,destination):
    """Copy every rational exponential coefficient; never evaluate exp."""
    return destination.EP({j:destination.ExpSum([(c,e) for e,c in v.terms.items()])
                           for j,v in value.c.items()})

class FixedProfile:
    """One chosen smooth rhohat with signed mass one; H_j >= sup|rhohat^(j)|/j!.

    The class computes bounds, not the profile or a check of its exact mass.
    No current field, covariance or signed target is needed to define rho.
    """
    def __init__(self,*,hull,Rmin,Rmax,width,profile):
        if set(hull)!={'l','u','Xa','Xb','t1'}:raise ValueError('Whole moving profile hull required')
        self.hull={k:exact(v) for k,v in hull.items()}
        if self.hull['t1']!=width:raise ValueError('Keep the original flat width')
        self.geometry=pipeline.fresh.geometry.geometry(**self.hull)
        self.Rmin,self.Rmax=exact(Rmin),exact(Rmax)
        if not 0<self.Rmin<self.Rmax:raise ValueError('Positive common norm hull required')
        if self.geometry['Rminus']<self.Rmin or self.geometry['Rplus']>self.Rmax:
            raise ValueError('Moving profile support must fit the whole norm hull')
        self.profile=tuple(exact(v) for v in profile)
        if not self.profile or min(self.profile)<0:raise ValueError('Nonnegative finite profile ledger required')
        self.max_order=len(self.profile)-1
        self.length=self.hull['u']-self.hull['l']
        if self.length*self.profile[0]<1:raise ValueError('Height/support cannot contain unit signed mass')
        self.cstar=self.length*self.hull['u']**2*self.profile[0]/2
        self.coordinates=signed.projector.target.ChartBudget(1,1)
        self.floor=EP(signed.maps.ExpSum([(1,self.geometry['log_weight_comparison'])]))
    def check(self,n):
        order(n)
        if n>self.max_order:raise ValueError('Outside the fixed profile ledger')
    @lru_cache(None)
    def dilation(self,n):
        self.check(n)
        if not n:return Q(2)
        return sum((Q(j,2*n)*self.coordinates.log_s(j)*self.dilation(n-j)
                    for j in range(1,n+1)),Q(0))
    def increment(self,j):
        self.check(j)
        return Q(0) if not j else self.Rmax*self.dilation(j)+self.dilation(j-1)
    @lru_cache(None)
    def power(self,n,k):
        self.check(n);order(k)
        if not k:return Q(int(n==0))
        return sum((self.increment(j)*self.power(n-j,k-1) for j in range(1,n+1)),Q(0)) if n>=k else Q(0)
    @lru_cache(None)
    def rho(self,n):
        self.check(n)
        return sum((self.dilation(i)*self.profile[k]*self.power(j,k)
                    for i in range(n+1) for j in range(n+1-i) for k in range(j+1)),Q(0))
    def weighted_rho(self,n):return self.floor*self.rho(n)
    def c_rho(self,n):
        self.check(n)
        return self.cstar*(4+sum(32**i*64**(j-i) for j in range(1,n+1) for i in range(j+1)))
    def sigma(self,n,e):
        self.check(n)
        if e not in (1,2):raise ValueError('Distinct axial/angular channel required')
        inverse=sum(Q(comb(e+j-1,j))*self.Rmin**(-e-j) for j in range(n+1))
        return self.length*inverse*sum(self.dilation(j) for j in range(n+1))*self.rho(n)**2

class Source(signed.projector.ProjectedSource):
    """Adapt the SAME Q35 norm objects; no synthetic covariance/current rows."""
    def __init__(self,haar,profile):
        if not isinstance(haar,pipeline.haar.InitialHaar):raise TypeError('Q35 initial prefix required')
        if getattr(haar.c_rho,'__self__',None) is not profile:raise ValueError('Whole pressure moment must use this fixed profile')
        primary=haar.initial.budget.primary
        if (primary.chart.Rmin,primary.chart.Rmax,primary.primary.root.width)!=(profile.Rmin,profile.Rmax,profile.hull['t1']):
            raise ValueError('Profile and original root/chart must agree')
        self.haar,self.profile=haar,profile
        self.geometry=profile.geometry
        self.pull=signed.projector.LogPullback(profile.Rmax)
    def bump(self,n,e):return self.profile.floor*self.pull.pullback(n,e)*self.profile.sigma(n,e)
    def source(self,n,e):
        f=pipeline.fresh
        out=self.haar.projected_source(n,e,
            pullback=lambda j,k:rebind(self.pull.pullback(j,k),f),
            bump=lambda j,k:rebind(self.bump(j,k),f))
        for key in ('row','moment','H'):out[key]=rebind(out[key],signed.projector.haar.sharp.fresh)
        return dict(out,fixed_rho_squared_choice=True,actual_fixed_profile_premises_required=True)

def build(root,stress,*,profile,profile_hull,**kwargs):
    """Q35 arguments except c_rho; fixed supplies all its old keys EXCEPT rho.

    The specified rho/profile must be the one used by the actual fixed maps.
    This computes fresh same-choice norm recipes, not an historical identity.
    """
    if 'c_rho' in kwargs:raise ValueError('Independent moving pressure callback is not allowed')
    m=kwargs['max_output'];order(m)
    p=FixedProfile(hull=profile_hull,Rmin=kwargs['Rmin'],Rmax=kwargs['Rmax'],width=root.width,profile=profile)
    if p.max_order<m+13:raise ValueError('Q35 fixed-data calls require profile jets through m+13')
    remaining=kwargs['fixed']
    def fixed(n,pole):
        d=remaining(n,pole)
        if set(d)!=set(pipeline.fresh.FIXED)-{'rho'}:raise ValueError('Supply precisely the remaining fixed-map inputs')
        return dict(d,rho=rebind(p.weighted_rho(n),pipeline.fresh))
    haar=pipeline.build(root,stress,**dict(kwargs,fixed=fixed),c_rho=p.c_rho)
    return signed.SignedTarget(Source(haar,p),m)

def record():
    return dict(result='conditional-fixed-profile-to-original-signed-target',
        independent_rho_c_rho_projector_callbacks_removed=True,
        same_choice_recomputed=True,historical_unspecified_projector_identity=False,
        profile_order='m+13',covariance_order='m+13',common_pole='12*(m+13)+18',
        native_target_pole='2*p+3*n',complete_curl_needs_next_target_jet=True,
        full_moving_pressure_and_two_edge_primitive=True,
        actual_profile_fixed_map_root_premises_required=True,
        full_signed_covariance_fresh_means_source_or_residual_certified=False,
        blowup=False)
