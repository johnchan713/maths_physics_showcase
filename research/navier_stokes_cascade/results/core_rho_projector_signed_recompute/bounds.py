"""Fresh same-choice signed-source/root/covariance/pressure norm composition."""
from fractions import Fraction as Q
from math import ceil,comb,factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('recomputed_rho_family',HERE.parent/'core_rho_squared_projector_realization/bounds.py')
family_module=importlib.util.module_from_spec(sp);sp.loader.exec_module(family_module)
moving=family_module.moving;real=moving.real;signed=real.signed;projector=signed.projector
EP,order=real.EP,real.order

def signature(family):
 c=family.compensator
 return (tuple(sorted(c.hull.items())),tuple(c.profile),c.Rmin,c.Rmax,c.update.target.width)

class ChosenSource(projector.ProjectedSource):
 def __init__(self,haar,family):
  self.family=family;self.fingerprint=signature(family)
  super().__init__(haar,family.native,dict(family.compensator.hull))
 def bump(self,n,e):
  if signature(self.family)!=self.fingerprint:raise ValueError('Chosen pressure-profile/projector data changed')
  return super().bump(n,e)
 def source(self,n,e):
  r=super().source(n,e)
  return dict(r,projector_definition='sigma_e=R^-e sqrt(s) rho^2/int rhohat^2',
   actual_fixed_bump_seminorms_still_required=False,
   original_fixed_rho_profile_mass_support_and_jets_still_required=True)

class RecomputedSigned:
 """Recompute conditional bounds; no actual continuum field evaluation.

 Initial Haar rows/fields are the SAME pre-signed prefix. This installs
 the explicit new sigma choice in its projected-source recipe, then
 constructs NEW signed-target/root/curl and complete covariance objects.
 """
 def __init__(self,haar,family,max_output):
  if not isinstance(family,family_module.RhoSquaredProjectors):raise TypeError('Explicit rho-squared projector definition required')
  if not isinstance(haar,projector.haar.InitialHaar):raise TypeError('Complete original pre-signed Haar recipe required')
  order(max_output);self.max_output=max_output;self.family=family;self.haar=haar
  c=family.compensator;self.fingerprint=signature(family)
  if haar.initial.budget.primary is not c.update.primary:raise ValueError('Projectors and initial rows must use the identical primary chart/root')
  self.primary=haar.initial.budget.primary;self.original_root=self.primary.primary.root
  if haar.max_output<max_output+3 or c.max_order<max_output+4:raise ValueError('Need source/target through n+3 and the whole initial pressure moment through n+4')
  self.source=ChosenSource(haar,family)
  self.target=signed.SignedTarget(self.source,max_output+3)
  self.update=real.RealCovariance(self.target,max_output)
  self.compensator=moving.MovingCompensator(self.update,dict(c.hull),tuple(c.profile))
  if signature(family_module.RhoSquaredProjectors(self.compensator))!=self.fingerprint:raise ValueError('The pressure update must retain the same chosen profile and hull')
 def check(self):
  if signature(self.family)!=self.fingerprint:raise ValueError('Chosen projector profile or geometry changed after recomputation')
  if signature(family_module.RhoSquaredProjectors(self.compensator))!=self.fingerprint or self.compensator.update is not self.update:raise ValueError('Recomputed pressure/profile identity changed')
  if self.target.source is not self.source or self.source.family is not self.family or self.update.target is not self.target:raise ValueError('Signed-source/root/covariance choice identity changed')
  if self.source.haar is not self.haar or getattr(self.source.sigma_norm,'__self__',None) is not self.family:raise ValueError('Chosen source or projector recipe replaced')
  if self.haar.initial.budget.primary is not self.primary or self.target.root is not self.original_root or self.primary.primary.root is not self.original_root:raise ValueError('The unchanged positive primary root was replaced')
 def kernel(self,p):
  lo,hi=self.target.Llo,self.target.Lhi;r=min(Q(1),lo/4);w=self.target.width
  a,b=w*w,Q(4);q=max(0,ceil((p-3)/2))
  def D(c):return sum((Q(comb(q,j)*factorial(j))*(r*r/c)**j for j in range(q+1)),Q(0))/(2*c)
  return EP(moving.maps.ExpSum([(D(a),b/(lo-r)**2),(D(b),a/(lo-r)**2),
    (hi*moving.maps.edge_absorption(p,w),(a+b)/(r*r))]))
 def maps(self,n,p):
  self.check();a,b=self.compensator.Rmin,self.compensator.Rmax
  kp=1+b*Q(1,2)*self.kernel(p)
  def target(e):
   pos=sum(Q(comb(e,j))*b**(e-j) for j in range(min(e,n)+1))
   neg=sum(Q(comb(e+j-1,j))*a**(-e-j) for j in range(n+1))
   return pos*neg*kp
  return dict(order=n,pole=p,Kp=kp,primitive=target(1)+target(2),
   bump=self.family.weighted(n,1)+self.family.weighted(n,2))
 def covariance(self,n):self.check();return self.update.covariance(n)
 def pressure(self,n):
  self.check();return self.update.pressure(n,self.compensator.operators(self.maps))
 def ledger(self):
  return dict(projector_choice='rho-squared',new_source_target_curl_covariance_and_pressure_norms_recomputed=True,
   old_signed_target_reused=False,source_profile_required_order=self.max_output+3,
   sufficient_fixed_profile_interface_order=self.max_output+4,
   actual_fixed_rho_and_initial_prefix_inputs_still_required=True,actual_continuum_fields_evaluated=False,
   fresh_fast_slow_or_complete_residual_certified=False)
