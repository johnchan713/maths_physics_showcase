"""Actual Haar R-to-log-chart jet bound and fixed projector norm routing."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
haar=load('log_pullback_haar','core_symbolic_initial_haar_join/bounds.py')
target=load('log_pullback_target','core_primary_target_cutoff_budget/bounds.py')
EP,order,exact=haar.EP,haar.order,haar.sharp.exact

class LogPullback:
 """Three-variable Haar coefficient norms: R,Z,T -> x,Z,T, no torus jets."""
 def __init__(self,Rmax):
  self.Rmax=exact(Rmax)
  if self.Rmax<=0:raise ValueError('Same positive whole radial norm hull required')
  self.coordinates=target.ChartBudget(1,1) # only universal log_s jets are used
 def log_s(self,j):return self.coordinates.log_s(j)
 @lru_cache(None)
 def exponential(self,b,n):
  b=exact(b);order(n)
  if b<0:raise ValueError('Nonnegative exponential majorant parameter required')
  if not n:return Q(1)
  return sum((b*Q(j,n)*(int(j==1)+self.log_s(j))*self.exponential(b,n-j) for j in range(1,n+1)),Q(0))
 def increment(self,j):
  order(j)
  return Q(0) if not j else 2*int(j==1)+self.Rmax*self.exponential(Q(1,2),j)
 @lru_cache(None)
 def power(self,n,k):
  order(n);order(k)
  if not k:return Q(int(n==0))
  if n<k:return Q(0)
  return sum((self.increment(j)*self.power(n-j,k-1) for j in range(1,n+1)),Q(0))
 def composition(self,n):
  order(n);return sum((self.power(n,k) for k in range(n+1)),Q(0))
 def pullback(self,n,e):
  order(n)
  if e not in (1,2):raise ValueError('Original axial/angular weight required')
  v=e+1
  return EP(self.Rmax**v*Q(1,2)*sum((self.exponential(Q(v,2),i)*self.composition(j)
       for i in range(n+1) for j in range(n+1-i)),Q(0)))

class ProjectedSource:
 """Discharge PB, route the ACTUAL fixed bump's native seminorms and support.

 sigma_norm(n,e) bounds through-order R,Z,T jets uniformly in the band;
 it is NOT inferred from mass, positivity or compact support alone.
 """
 def __init__(self,initial_haar,sigma_norm,hull):
  if not isinstance(initial_haar,haar.InitialHaar):raise TypeError('Same-family complete initial Haar join required')
  if not callable(sigma_norm):raise TypeError('Original projector fixed-family seminorm recipe required')
  if set(hull)!={'l','u','Xa','Xb','t1'}:raise ValueError('Whole actual projector support hull required')
  self.haar,self.sigma_norm=initial_haar,sigma_norm
  primary=initial_haar.initial.budget.primary
  if exact(hull['t1'])!=primary.primary.root.width:raise ValueError('Retain the same original flat width')
  self.pull=LogPullback(primary.chart.Rmax)
  fresh=haar.sharp.fresh
  self.geometry=fresh.geometry.geometry(**hull)
  if self.geometry['Rminus']<primary.chart.Rmin or self.geometry['Rplus']>primary.chart.Rmax:raise ValueError('Projector moving support must remain in the original norm hull')
  self.floor=EP(fresh.ExpSum([(1,self.geometry['log_weight_comparison'])]))
 def bump(self,n,e):return self.floor*self.pull.pullback(n,e)*EP(self.sigma_norm(n,e))
 def source(self,n,e):
  r=self.haar.projected_source(n,e,pullback=self.pull.pullback,bump=self.bump)
  return dict(r,logarithmic_pullback_derived=True,actual_fixed_bump_seminorms_still_required=True)
