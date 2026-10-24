"""Fixed-profile moving compensator and all-new-pole radial moment bounds."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
real=load('moving_pressure_real','core_symbolic_real_covariance_join/bounds.py')
EP,order,exact=real.EP,real.order,real.signed.exact
maps=real.signed.maps
geometry=real.signed.projector.haar.sharp.fresh.geometry

class MovingCompensator:
 """Same rho=s^-1/2 rhohat(R/sqrt(s)); fixed rhohat jets remain inputs.

 profile[k] bounds sup |rhohat^(k)|/k! on the entire fixed profile.
 Unit signed mass and actual support identification are exact premises.
 """
 def __init__(self,update,hull,profile):
  if not isinstance(update,real.RealCovariance):raise TypeError('Same-family real covariance update required')
  if set(hull)!={'l','u','Xa','Xb','t1'}:raise ValueError('Whole actual fixed compensator hull required')
  self.update=update;self.hull={k:exact(v) for k,v in hull.items()}
  if self.hull['t1']!=update.target.width:raise ValueError('Keep the original flat width')
  self.geom=geometry.geometry(**self.hull)
  chart=update.primary.chart;self.Rmin,self.Rmax=chart.Rmin,chart.Rmax
  if self.geom['Rminus']<self.Rmin or self.geom['Rplus']>self.Rmax:raise ValueError('Moving compensator must fit the same whole radial norm hull')
  self.profile=tuple(exact(v) for v in profile)
  if not self.profile or min(self.profile)<0:raise ValueError('Nonnegative fixed-profile jet bounds required')
  self.max_order=len(self.profile)-1
  if self.max_order<update.max_output+1:raise ValueError('Whole pressure moment requires the next compensator jet')
  l,u=self.hull['l'],self.hull['u']
  if (u-l)*self.profile[0]<1:raise ValueError('Declared height/support cannot contain a unit-mass compensator')
  self.cstar=(u-l)*u*u*self.profile[0]/2
  self.floor=EP(maps.ExpSum([(1,self.geom['log_weight_comparison'])]))
  self.chart=real.signed.projector.target.ChartBudget(1,1)
 def check(self,n):
  order(n)
  if n>self.max_order:raise ValueError('Outside fixed-profile jet ledger')
 @lru_cache(None)
 def dilation(self,n):
  """Total-degree coefficient majorant for s^-1/2, with no x-log jet."""
  self.check(n)
  if not n:return Q(2)
  return sum((Q(j,2*n)*self.chart.log_s(j)*self.dilation(n-j) for j in range(1,n+1)),Q(0))
 def increment(self,j):
  self.check(j)
  return Q(0) if not j else self.Rmax*self.dilation(j)+self.dilation(j-1)
 @lru_cache(None)
 def power(self,n,k):
  self.check(n);order(k)
  if not k:return Q(int(n==0))
  if n<k:return Q(0)
  return sum((self.increment(j)*self.power(n-j,k-1) for j in range(1,n+1)),Q(0))
 def rho_unweighted(self,n):
  self.check(n)
  return sum((self.dilation(i)*self.profile[k]*self.power(j,k)
      for i in range(n+1) for j in range(n+1-i) for k in range(j+1)),Q(0))
 def rho(self,n):return self.floor*self.rho_unweighted(n)
 def c_rho(self,n):
  self.check(n)
  return self.cstar*(4+sum(32**i*64**(j-i) for j in range(1,n+1) for i in range(j+1)))
 def moment(self,pole):
  """All radial weights -1,0,1,2, normalized Haar, fixed union hull."""
  pole=exact(pole)
  if pole<0:raise ValueError('Nonnegative source pole required')
  return (self.Rmax-self.Rmin)*(1/self.Rmin+1+self.Rmax+self.Rmax**2)*maps.edge_absorption(pole,self.hull['t1'])
 def operators(self,remaining):
  """Only Kp, primitive and bump remain new-domain lower hypotheses."""
  if not callable(remaining):raise TypeError('Remaining actual operator bounds required')
  def recipe(n,p):
   self.check(n);r=remaining(n,p)
   if set(r)!={'order','pole','Kp','primitive','bump'}:raise ValueError('Keep precisely the remaining new-domain map inputs')
   order(r['order']);order(r['pole'])
   if r['order']<n or r['pole']<p:raise ValueError('Remaining maps do not cover the new order/pole')
   return dict(r,rho=self.rho(n),Kmom=self.moment(p),c_rho=self.c_rho(n))
  return recipe
 def pressure(self,n,remaining):return self.update.pressure(n,self.operators(remaining))
 def ledger(self):
  return dict(moving_compensator_and_radial_moment_recipes_derived=True,
   actual_fixed_profile_support_mass_and_seminorms_required=True,
   pressure_primitive_projector_and_actual_sigma_family_still_required=True,
   full_residual_or_blowup_certified=False)
