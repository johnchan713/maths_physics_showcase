"""Current Q38 scalars plus complete fast changes, then the same slow response."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
maps=load('q39_maps','core_current_head_q38_new_pole_maps/new_maps.py')
post=load('q39_post','core_post_signed_mean_closure/budgets.py')
EP,order=maps.EP,maps.order
fresh=maps.update.signed.projector.haar.sharp.fresh

def convert(v):
 if hasattr(v,'c') and any(hasattr(c,'terms') for c in v.c.values()):
  return maps.update.profile.rebind(v,fresh)
 return EP(v)

class CurrentScalars:
 """Actual new fast norm recipe remains conditional; no scalar/old-mean callback.

 fast(n) supplies through-n weighted norms for complete x,u,z, order and
 pole. Exact zero Haar and the original fast-map identity remain premises.
 """
 def __init__(self,current,fast,max_output):
  if not isinstance(current,maps.NewPoleMaps):raise TypeError('Same Q38 pressure/profile state required')
  order(max_output)
  if current.state.max_output<max_output+1:raise ValueError('Slow radial response needs the next current scalar')
  if not callable(fast):raise TypeError('Actual complete fast-field norm recipe required')
  self.current,self.fast,self.max_output=current,fast,max_output
  self.initial=current.target.source.haar.initial
  cp=current.profile
  self.moment=(cp.Rmax-cp.Rmin)*(1/cp.Rmin+1+cp.Rmax+cp.Rmax**2)
 def check(self,n):
  order(n)
  if n>self.max_output+1:raise ValueError('Outside current scalar ledger')
 @lru_cache(None)
 def old_full(self,n):
  order(n);self.check(max(0,n-1))
  b=self.initial.budget;f=b._fast(n);s=self.initial.slow_fields(n,False)
  edge=maps.update.signed.maps.edge_absorption(b.pole,self.current.target.width)
  return {k:edge*convert(f[k])+convert(s[j]) for k,j in [('x','X'),('u','U'),('z','Z')]}
 @lru_cache(None)
 def fast_unweighted(self,n):
  order(n);self.check(max(0,n-1));v=self.fast(n)
  if set(v)!={'order','pole','x','u','z'}:raise ValueError('Complete fast fields and finite domain declaration required')
  order(v['order']);order(v['pole'])
  if v['order']<n:raise ValueError('Fast field ledger too short')
  edge=maps.update.signed.maps.edge_absorption(v['pole'],self.current.target.width)
  return {k:edge*convert(v[k]) for k in ('x','u','z')}
 @lru_cache(None)
 def postfast(self,n):
  self.check(n)
  current=self.current.pressure(n)['updated_scalars']
  change=post.scalar.fast_change(n,self.old_full(n+1),self.fast_unweighted(n+1),self.moment)
  return {k:current[k]+change[k] for k in ('P','Jtheta','Jz')}
 def slow(self,n,weighted=False):
  order(n)
  if n>self.max_output:raise ValueError('Outside slow response ledger')
  if type(weighted) is not bool:raise TypeError('Choose weighted or unweighted response')
  b=self.initial.budget
  k={key:convert(v) for key,v in b.slow(n).items()}
  kn={key:convert(v) for key,v in b.slow(n+1).items()}
  keys={'H','Kp','Kz','Kg','Kpsi'}
  if set(k)!=keys or set(kn)!=keys:raise ValueError('Same complete five-map response required')
  out=post.preserved.response(n,self.postfast(n),self.postfast(n+1),k,kn)
  if weighted:out={key:convert(b.floor)*v for key,v in out.items()}
  return out

def record():
 k=Q(1,100000);B=1-k;alpha=2-3*k;D=2-4*k
 return dict(result='current-scalar-fast-change-and-slow-response-join',
  old_full_means_derived=True,signed_scalar_change_added_once=True,
  complete_axial_fast_and_old_correlations_required=True,
  old_and_new_poles_unweighted_separately=True,
  postfast_order=str(min(Q(2),alpha,B+D,2*D)),fast_change_order=str(min(B+D,2*D)),
  Q38_pressure_ceiling='m+1',new_fast_jet_ceiling='m+2',Q36_target_ceiling='m+4',
  Q36_profile_covariance_ceiling='m+17',actual_fast_source_bound_required=True,
  postslow_reset_full_source_residual_or_blowup_certified=False)
