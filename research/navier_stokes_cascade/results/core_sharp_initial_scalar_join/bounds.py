"""Sharper INITIAL scalar identities joined to the actual symbolic norm route."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
fresh=load('sharp_join_fresh','core_symbolic_fresh_mean_join/bounds.py')
initial=load('sharp_join_initial','core_initial_scalar_seed_budget/proof_checks.py')
EP,order,exact=fresh.EP,fresh.order,fresh.exact

class SharpInitial:
 """No new scalar/current-row inputs. Original zero-Haar/five-map premises apply.

 fast_hull must contain the WHOLE primary, fast potential and A1 supports,
 with all input norms valid there. Reserved slow hull is separate.
 """
 def __init__(self,budget,max_output,fast_hull):
  if not isinstance(budget,fresh.FreshMeans):raise TypeError('Complete same-family primary/fast symbolic recipe required')
  order(max_output)
  if budget.max_output<max_output+1:raise ValueError('Reserve the next scalar output for the whole moving axial moment')
  a,d=map(exact,fast_hull)
  if not 0<a<d:raise ValueError('Whole positive fast potential/cutoff hull required')
  if a<budget.primary.chart.Rmin or d>budget.primary.chart.Rmax:raise ValueError('Fast support hull must lie in the original norm hull')
  self.budget,self.max_output=budget,max_output
  self.M=initial.coupled.geometry(a,d)
  self.Ms=initial.coupled.geometry(budget.hull['Rminus'],budget.hull['Rplus'])
  width=budget.primary.primary.root.width
  self.edge=fresh.maps.edge_absorption(budget.pole,width)
 @lru_cache(None)
 def initial_scalars(self,n):
  self._scalar_order(n);M=self.M
  w=EP(self.budget.primary.covariance_R(n)['W'])*self.edge
  wn=EP(self.budget.primary.covariance_R(n+1)['W'])*self.edge
  return dict(P=2*M[-1]*w+(n+1)*M[0]*wn,Jtheta=M[2]*w,
              Jz=2*M[1]*w+Q(n+1,2)*M[2]*wn)
 @lru_cache(None)
 def fast_scalar_change(self,n):
  self._scalar_order(n);M=self.M;f=self.budget._fast(n+1)
  T=self.edge*f['u'];X=self.edge*f['x']
  # Same historical zero-old-Haar fast quadratic recipe. 2mu>=1.
  out=initial.coupled.fast_budget(n,M,0,0,T,X)
  return dict(P=out['P'],Jz=out['Jz'],Jtheta=M[2]*T*T)
 @lru_cache(None)
 def preslow(self,n):
  self._scalar_order(n);a=self.initial_scalars(n);b=self.fast_scalar_change(n)
  return {k:a[k]+b[k] for k in a}
 @lru_cache(None)
 def slow_fields(self,n,weighted=False):
  order(n)
  if n>self.max_output+3:raise ValueError('Outside the declared slow field ledger')
  if type(weighted) is not bool:raise TypeError('Choose the original unweighted or actual-hull weighted field norm')
  a,b=self.preslow(n),self.preslow(n+1)
  k,kn=self.budget.slow(n),self.budget.slow(n+1)
  if any(set(v)!=set(('H','Kp','Kz','Kg','Kpsi')) for v in (k,kn)):raise ValueError('Every original slow response norm required')
  k={key:EP(v) for key,v in k.items()};kn={key:EP(v) for key,v in kn.items()}
  v=dict(U=k['H']*(k['Kp']*a['P']+k['Kz']*a['Jz']),Z=k['H']*k['Kg']*a['Jtheta'],
         X=(n+1)*kn['H']*kn['Kpsi']*b['Jtheta'])
  return {key:self.budget.floor*x for key,x in v.items()} if weighted else v
 @lru_cache(None)
 def output(self,n):
  order(n)
  if n>self.max_output+1:raise ValueError('Outside the declared finite scalar/moment ledger')
  v=self.slow_fields(n+2);b=EP(self.budget.fixed(n+2,2*self.budget.pole)['b'])
  # Existing C=1 initial five-map remainder, G=0 on the reserved hull.
  return initial.coupled.slow_budget(n,self.Ms,b,0,0,v['X'],v['U'],v['Z'])
 def axial_moment(self,n,c_rho):
  order(n)
  if n>self.max_output:raise ValueError('Outside the declared finite scalar/moment ledger')
  v=self.output(n+1)
  return (n+1)*(v['Jz']+EP(c_rho)*v['P'])
 def ledger(self):
  return dict(preslow_scalar_order=1,slow_Haar_tangential_order=1,slow_Haar_radial_order=2,
              postslow_scalar_order=2,whole_axial_moment_order=3,
              covariance_order_for_scalar=self.max_output+12,covariance_order_for_axial=self.max_output+13,
              initial_zero_Haar_and_reserved_five_map_premises_required=True,
              later_complex_prefix_or_full_residual_certified=False)
 def _scalar_order(self,n):
  order(n)
  if n>self.max_output+4:raise ValueError('Outside the declared fresh scalar input ledger')
