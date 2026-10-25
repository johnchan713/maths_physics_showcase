"""Complete first-real postfast pressure/scalar refresh and sharp slow input.

All coefficient bounds remain conditional on the inherited actual-field,
profile, support and operator contracts.  This module does not turn the
finite positive-polynomial program into a continuum enclosure.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

previous=load('postfast_previous','core_initial_quotient_fast_inverse/bounds.py')
refresh=load('postfast_refresh','core_real_refresh_obligation_join/bounds.py')
potential=load('postfast_potential','core_axial_potential_pressure_gain/bounds.py')
EP,order,exact=previous.EP,previous.order,previous.exact
SCALARS=('P','Jtheta','Jz')
FIELDS=('x','u','z')
SLOW=('H','Kp','Kz','Kg','Kpsi')

def complete_pressure_envelope(prefast,initial_fast,initial_slow,signed_wave,new_fast):
 """Add the five sequential pressure states, with no pressure reset."""
 parts=dict(prefast=EP(prefast),initial_fast=EP(initial_fast),
  initial_slow=EP(initial_slow),signed_wave=EP(signed_wave),new_fast=EP(new_fast))
 return dict(parts=parts,pm=sum(parts.values(),EP()),
  complete_pressure_memory_retained=True,additive_pressure_gauge_separate=True)

def postfast_scalar_refresh(n,pressure_result,old_full,fast_increment,moment):
 """Add only the new fast change to the already post-signed scalars."""
 order(n)
 if set(pressure_result['updated_scalars'])!=set(SCALARS):
  raise ValueError('All current post-signed scalar defects required')
 if set(old_full)!=set(FIELDS) or set(fast_increment)!=set(FIELDS):
  raise ValueError('Three complete old and new-fast mean components required')
 old={k:EP(v) for k,v in old_full.items()};inc={k:EP(v) for k,v in fast_increment.items()};M=EP(moment)
 x,u,z=(inc[k] for k in FIELDS);ox,ou,oz=(old[k] for k in FIELDS)
 S=2*ox*x+x*x;H=ox*z+oz*x+x*z;T=2*ou*u+u*u
 Zf=2*oz*z+z*z;K=oz*u+ou*z+u*z;q=M*(S+T+(n+1)*H)
 change=dict(P=q,Jtheta=M*K,Jz=M*Zf+q*Q(1,2))
 return {k:EP(pressure_result['updated_scalars'][k])+change[k] for k in SCALARS}

def sharp_slow_fields(n,current,current_next,slow,slow_next,floor=1):
 """The original weighted five-map evaluated on the refreshed scalars."""
 order(n)
 if set(current)!=set(SCALARS) or set(current_next)!=set(SCALARS):
  raise ValueError('Complete current scalar defects at n and n+1 required')
 if set(slow)!=set(SLOW) or set(slow_next)!=set(SLOW):
  raise ValueError('Complete original slow five-map at n and n+1 required')
 a={k:EP(v) for k,v in current.items()};an={k:EP(v) for k,v in current_next.items()}
 k={key:EP(v) for key,v in slow.items()};kn={key:EP(v) for key,v in slow_next.items()};w=EP(floor)
 return dict(u=w*k['H']*(k['Kp']*a['P']+k['Kz']*a['Jz']),
  z=w*k['H']*k['Kg']*a['Jtheta'],
  x=w*(n+1)*kn['H']*kn['Kpsi']*an['Jtheta'])

def a1_time_error_request(n,kappa=Q(1,100000)):
 """Strict current-potential request for the new exact A1 time error.

 A1 Q Ez=A1 Dz Q Phi has source order two.  The target is the first-real
 final-row order nu, so equality at p-1 is insufficient.
 """
 order(n);k=exact(kappa)
 if not 0<k<=Q(1,100):raise ValueError('Retained positive small-kappa regime required')
 source=Q(2);target=Q(3)-6*k
 p=max(0,(target-source)//k+1);J=n+1 if p==0 else n+2*p+4
 return dict(output_order=n,source_stage='post-signed/pre-new-fast',
  identity='A1 Qcirc Ez_postsigned=A1 Dz Qcirc Phi_postsigned',
  Phi='-2 G_R Psi+gamma^2+Wzz+pm',postfast_or_postslow_pressure_substituted=False,
  source_order=source,target_order=target,
  integration_order=p,raw_Phi_order=J,strict_margin=source+p*k-target,
  current_fields_for_potential=J+6,current_Ez_for_potential=J+4,
  coarse_pressure_for_potential=J+5,
  conservative_pressure_field_recipe=J+8,
  current_Etheta_for_fine_pressure=J+4,fixed_G_for_GR=J+1,
  previous_full_row_integration_order=99999,
  previous_full_row_raw_source_order=n+200001,
  direct_order_two_alone_reaches_target=source>target,
  actual_current_Phi_high_jets_and_support_materialized=False)

def rates(kappa=Q(1,100000)):
 k=exact(kappa)
 if not 0<k<=Q(1,100):raise ValueError('Retained positive small-kappa regime required')
 B=1-k;alpha=2-3*k;D=2-4*k;nu=3-6*k
 return dict(B=B,alpha=alpha,D=D,nu=nu,
  new_fast_tangential=D,new_fast_radial=D+1,
  new_fast_scalar_change=min(B+D,2*D),postfast_scalar=alpha,
  sharp_slow_tangential=alpha,sharp_slow_axial=alpha,sharp_slow_radial=alpha+1,
  postslow_scalar=alpha+1,whole_postslow_axial=alpha+2,
  final_fast_difference=D+1-2*k,
  final_fast_difference_strict_margin=D+1-2*k-nu,
  current_A1_source=Q(2),actual_current_A1_high_jets_materialized=False,
  complete_residual_improvement_or_blowup_certified=False)

class InitialPostfastRefresh:
 """One-way conditional join from the corrected fast inverse to fresh means."""
 def __init__(self,inverse):
  if not isinstance(inverse,previous.InitialQuotientFastInverse):
   raise TypeError('Corrected initial quotient fast inverse required')
  self.inverse=inverse;self.current=inverse.current;self.state=inverse.state
  self.budget=inverse.budget;self.pressure=inverse.pressure
 def _check(self):
  self.inverse._check()
  if self.current.state is not self.state or self.current.budget is not self.budget:
   raise ValueError('Inherited current state or original mean recipe changed')
 def _fixed(self,n):
  order(n);self._check();N=n+3
  base=self.budget.fixed(N,self.pressure.pole)
  if not {'b','V','G','c'}<=set(base):raise ValueError('Original fixed coefficients required')
  chart=self.state.primary.chart;lo=chart.Rmin
  out={k:EP(base[k]) for k in ('b','V','G','c')}
  out.update(graph=EP(2+6*chart.B*sum((j+1)*lo**(-j) for j in range(N+1))),
   r=EP(sum(lo**(-j-1) for j in range(N+1))),rho=self.state.compensator.rho(n+1),
   Kp=self.pressure.KT0,Kmom=EP(self.state.compensator.moment(self.pressure.pole)))
  return out
 def old_full(self,n):
  order(n);self._check();f=self.budget._fast(n);s=self.current.initial.slow_fields(n,True)
  return dict(x=EP(f['x'])+EP(s['X']),u=EP(f['u'])+EP(s['U']),z=EP(f['z'])+EP(s['Z']))
 def new_fast(self,n):
  f=self.inverse.fast_fields(n);return {k:EP(f[k]) for k in FIELDS}
 def prefast_pressure(self,n):
  order(n);fixed=self._fixed(n);N=n+3;W=EP(self.state.primary.covariance_R(N)['W']);d=n+4
  gr=(d*fixed['graph']+d+2*fixed['r'])*W;P=fixed['Kmom']*gr
  return dict(gr=gr,P=P,pm=fixed['Kp']*(gr+fixed['rho']*P),input_order=N)
 def fast_difference(self,n):
  order(n);return previous.prior.fast.difference(n,self.old_full(n+3),self.new_fast(n+3),self._fixed(n))
 def pressure_envelope(self,n):
  order(n);old=self.current.rows(n);base=self.prefast_pressure(n);new=self.fast_difference(n)
  out=complete_pressure_envelope(base['pm'],old['original_fast_difference']['pm'],
   old['original_slow_difference']['pm'],old['complete_signed_covariance_difference']['pm'],new['pm'])
  out.update(prefast=base,new_fast_difference=new,input_order=n+3)
  return out
 def scalars(self,n):
  order(n);self._check();result=self.state.pressure(n);N=n+1
  return postfast_scalar_refresh(n,result,self.old_full(N),self.new_fast(N),
   EP(self.state.compensator.moment(self.pressure.pole)))
 def slow_fields(self,n):
  order(n);self._check()
  return sharp_slow_fields(n,self.scalars(n),self.scalars(n+1),
   self.budget.slow(n),self.budget.slow(n+1),self.budget.floor)
 def recipe(self,n):
  order(n);q=a1_time_error_request(n)
  return dict(output_order=n,pressure_components=('prefast','initial_fast','initial_slow','signed_wave','new_fast'),
   scalar_recipe='post-signed updated scalars plus only the new actual-fast covariance change',
   slow_recipe='same original moving-hull weighted five-map on refreshed scalars',
   new_A1=q,next_obligation='materialize the current Phi high jets/support and exact A1 term, then assemble the complete final row including the sharp slow difference')
 def rates(self):return rates()
