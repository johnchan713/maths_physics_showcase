"""Native Haar pressure/target primitives at the complete new covariance pole."""
from fractions import Fraction as Q
from functools import lru_cache
from math import ceil,comb,factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('native_haar_moving',HERE.parent/'core_moving_pressure_compensator_join/bounds.py')
moving=importlib.util.module_from_spec(sp);sp.loader.exec_module(moving)
EP,order,exact=moving.EP,moving.order,moving.exact

class NativeHaarPrimitives:
 def __init__(self,compensator):
  if not isinstance(compensator,moving.MovingCompensator):raise TypeError('Same moving compensator and covariance state required')
  self.compensator=compensator;self.target=compensator.update.target
  if compensator.hull['Xb']/compensator.hull['Xa']!=self.target.ratio:raise ValueError('Retain the same annular log length')
  self.a,self.b=compensator.Rmin,compensator.Rmax
 @lru_cache(None)
 def kernel(self,p):
  p=exact(p)
  if p<0:raise ValueError('Nonnegative finite source pole required')
  lo,hi=self.target.Llo,self.target.Lhi;r=min(Q(1),lo/4)
  w=self.target.width;a,b=w*w,Q(4);q=max(0,ceil((p-3)/2))
  def D(c):return sum((Q(comb(q,j)*factorial(j))*(r*r/c)**j for j in range(q+1)),Q(0))/(2*c)
  return EP(moving.maps.ExpSum([(D(a),b/(lo-r)**2),(D(b),a/(lo-r)**2),
      (hi*moving.maps.edge_absorption(p,w),(a+b)/(r*r))]))
 def pressure_constant(self,n,p):
  self.compensator.check(n)
  return EP(1)+self.b*Q(1,2)*self.kernel(p)
 def target_constant(self,n,p,e):
  self.compensator.check(n)
  if e not in (1,2):raise ValueError('Original axial/angular channel required')
  positive=sum(Q(comb(e,j))*self.b**(e-j) for j in range(min(e,n)+1))
  negative=sum(Q(comb(e+j-1,j))*self.a**(-e-j) for j in range(n+1))
  return positive*negative*self.pressure_constant(n,p)
 def operators(self,bump):
  if not callable(bump):raise TypeError('Actual fixed sigma-family bump bound still required')
  def remaining(n,p):
   self.compensator.check(n)
   return dict(order=n,pole=p,Kp=self.pressure_constant(n,p),
    primitive=self.target_constant(n,p,1)+self.target_constant(n,p,2),bump=EP(bump(n,p)))
  return self.compensator.operators(remaining)
 def pressure(self,n,bump):return self.compensator.update.pressure(n,self.operators(bump))
 def ledger(self):
  return dict(Haar_only=True,native_RZT_derivative_loss=0,new_pole_pressure_and_target_primitives_derived=True,
   exact_pointwise_in_slow_variables_zero_mass_required=True,actual_sigma_and_fixed_profile_support_inputs_required=True,
   oscillatory_compact_maps_or_full_residual_certified=False)
