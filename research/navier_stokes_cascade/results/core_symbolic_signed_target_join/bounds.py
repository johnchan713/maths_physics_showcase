"""Two-edge primitive to the original signed root, with symbolic constants."""
from fractions import Fraction as Q
from functools import lru_cache
from math import ceil,comb,factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('symbolic_signed_projector',HERE.parent/'core_log_pullback_projector_join/bounds.py')
projector=importlib.util.module_from_spec(sp);sp.loader.exec_module(projector)
EP,order,exact=projector.EP,projector.order,projector.exact
maps=projector.haar.sharp.fresh.maps

class SignedTarget:
 """Conditional norm composition only; actual bump-family validity remains open."""
 def __init__(self,source,max_order):
  if not isinstance(source,projector.ProjectedSource):raise TypeError('Same-family projected initial source required')
  order(max_order)
  if max_order>source.haar.max_output:raise ValueError('Projected-source derivative ledger too short')
  self.source,self.max_order=source,max_order
  primary=source.haar.initial.budget.primary
  self.root=primary.primary.root
  self.width=self.root.width;self.pole=2*source.haar.pole
  self.Rmin=primary.chart.Rmin;self.Ufloor=min(Q(1),self.Rmin**3)
  # Recover the SAME shell ratio from the recorded geometric gap identities.
  g=source.geometry;l=2*g['Rminus'];u=g['Rplus']/2
  da,db=g['left_log_gap_lower'],g['right_log_gap_lower']
  ra,rb=(2+da)/(2-da),(2+db)/(2-db)
  self.ratio=ra*rb*(u/l)**2
  self.Llo=2*(self.ratio-1)/(self.ratio+1);self.Lhi=self.ratio-1
  r=min(Q(1),self.Llo/4);a,b=self.width**2,Q(4);q=max(0,ceil((self.pole-3)/2))
  def D(c):return sum((Q(comb(q,j)*factorial(j))*(r*r/c)**j for j in range(q+1)),Q(0))/(2*c)
  self.primitive=EP(maps.ExpSum([(D(a),b/(self.Llo-r)**2),(D(b),a/(self.Llo-r)**2),
       (self.Lhi*maps.edge_absorption(self.pole,self.width),(a+b)/(r*r))]))
 def check(self,n):
  order(n)
  if n>self.max_order:raise ValueError('Outside the declared signed-target jet ledger')
 @lru_cache(None)
 def H(self,n,e):self.check(n);return self.source.source(n,e)['H']
 @lru_cache(None)
 def F(self,n,e):
  self.check(n);return self.primitive*self.H(n,e)+(self.H(n-1,e) if n else EP())
 @lru_cache(None)
 def K(self,n,e):
  self.check(n)
  return sum((self.Rmin**(-e)*self.source.pull.exponential(Q(e,2),j)*self.F(n-j,e)
              for j in range(n+1)),EP())
 @lru_cache(None)
 def reciprocal_weight(self,n):
  order(n)
  if not n:return Q(1)
  return sum((Q(j,n)*(self.width**2+4)*(j+1)*self.reciprocal_weight(n-j) for j in range(1,n+1)),Q(0))
 @lru_cache(None)
 def Vx(self,n,e):
  self.check(n)
  return sum((self.K(j,e)*self.reciprocal_weight(n-j) for j in range(n+1)),EP())
 def increment(self,j):
  order(j)
  return Q(0) if not j else 2*int(j==1)+Q(6,j)/self.Ufloor**j+self.source.pull.log_s(j)
 @lru_cache(None)
 def power(self,n,k):
  order(n);order(k)
  if not k:return Q(int(n==0))
  if n<k:return Q(0)
  return sum((self.increment(j)*self.power(n-j,k-1) for j in range(1,n+1)),Q(0))
 @lru_cache(None)
 def native(self,n,e):
  self.check(n)
  if e not in (1,2):raise ValueError('Original axial or angular signed channel required')
  return sum((self.Vx(k,e)*self.power(n,k) for k in range(n+1)),EP())
 def total(self,n):return self.native(n,1)+self.native(n,2)
 @lru_cache(None)
 def signed_amplitude(self,n):
  self.check(n)
  def V(j):return sum((EP(self.root.K(k))*self.total(j-k) for k in range(j+1)),EP())
  return self.root.mass*Q(1,2)*EP({1:1})*sum((EP(self.root.E(i))*V(j)*EP(self.root.reciprocal(k))*EP(self.root.data('wave',n-i-j-k))
      for i in range(n+1) for j in range(n-i+1) for k in range(n-i-j+1)),EP())
 def signed_curl(self,n):
  self.check(n+1)
  return sum((EP(self.root.data('curl',j))*(self.signed_amplitude(n-j)+(n-j+1)*self.signed_amplitude(n-j+1)) for j in range(n+1)),EP())
 def ledger(self):
  return dict(native_input_pole_offset=self.pole,native_pole_slope=3,
   signed_principal_order=Q(3,2)-Q(3,100000),signed_curl_order=2-Q(4,100000),
   actual_projector_uniform_seminorms_still_required=True,original_root_and_value_premises_still_required=True,
   actual_signed_wave_full_covariance_or_residual_certified=False)
