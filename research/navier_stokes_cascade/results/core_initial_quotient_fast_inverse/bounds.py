"""Correct the original prefast A1 tail and feed the whole source to the fast map.

Every estimate remains conditional on the original weighted field, profile,
support, pressure and moment inputs.  No continuum enclosure is manufactured.
"""
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

prior=load('initial_quotient_prior','core_initial_cancellation_fast_source/bounds.py')
EP,order,exact=prior.EP,prior.order,prior.exact

def prefast_quotient_rows(n,W,reciprocal,KT0,KA1,Kmom=0,rho=0):
 """Direct normalized A1 route for the ORIGINAL zero-mean prefast row.

 With S=Wrr, H=Wrz and C=Wtt,

   gr=r0+D_r(-S),  r0=-S/R-Dz H+C/R,
   pm=-S+T0(r0-rho P),  P=int Haar(r0).

 Hence A1 Q Ez=A1 Dz Q[Wzz-S+T0(r0-rho P)].  Q(rho P)=0,
 but the full P and pressure memory are returned to keep that distinction
 explicit.  W is a common normalized bound for all six covariance entries.
 """
 order(n);W,r,kp,ka,km,rr=map(EP,(W,reciprocal,KT0,KA1,Kmom,rho));d=n+4
 r0=(d+2*r)*W
 P=km*r0
 full_nonlocal=kp*(r0+rr*P)
 # Q costs two.  The compensator disappears only in this Q pressure channel.
 q_pressure=2*kp*r0
 q_phi=4*W+q_pressure
 return dict(r0=r0,P=P,full_local_pressure=W,full_nonlocal_pressure=full_nonlocal,
  Q_pressure=q_pressure,Q_Phi=q_phi,A1=d*ka*q_phi,
  covariance_input_order=n+3,pressure_order=n+1,operator_order=n,
  explicit_epsilon_power=2,
  exact_identity='A1 Qcirc Ez_prefast=A1 Dz Qcirc[Wzz-Wrr+T0(r0-rho P)]',
  complete_derivative_support_contract_required=True,
  actual_profile_support_or_field_enclosure_certified=False)

def fast_majorants(m,Etheta,Ez,Ez_radial,compact,compact_next):
 """Original weighted temporal/compact inverse on an ALREADY Q source.

 The source coefficients already contain Q's norm, so there is no second
 projection factor.  The output is the actual (x,u,z)=(beta,v,gamma) recipe,
 including gamma=gd-A1 gd and the external exact S factor.
 """
 order(m)
 if set(compact)!=set(('KT','KA','I')) or set(compact_next)!=set(('KT','KA','I')):
  raise ValueError('Complete original KT1, KA1 and time-floor recipes required')
 et,ez,ezr=map(EP,(Etheta,Ez,Ez_radial));op={k:EP(v) for k,v in compact.items()};nxt={k:EP(v) for k,v in compact_next.items()}
 S=EP({1:1});gain=comb(m+5,5)*24*2**(m+4)*op['I']*S
 tang=factorial(m+4)*gain*(et+(1+op['KA'])*ez)
 radial=factorial(m+5)*comb(m+5,5)*24*2**(m+5)*op['I']*S*nxt['KT']*ezr
 return dict(x=radial,u=tang,z=tang,actual_axial_velocity='gd-A1 gd',
  source_already_Q_projected=True,external_S_power=1,
  tangential_source_order=m+4,radial_source_order=m+5)

class InitialQuotientFastInverse:
 """Corrective one-way join from the prior current source to its fast map."""
 def __init__(self,current):
  if not isinstance(current,prior.CurrentFastSource):raise TypeError('Previous exact initial-cancellation source required')
  if not callable(getattr(current.budget,'compact',None)):raise TypeError('Original compact KT1/KA1/time-floor recipe required')
  self.current=current;self.state=current.state;self.budget=current.budget
  self.pressure=current.pressure;self.initial_pole=current.initial_pole
 def _check(self):
  if self.current.state is not self.state or self.current.budget is not self.budget or self.budget.pole!=self.initial_pole:
   raise ValueError('Original initial prefix, state or proof pole changed')
  self.state.check()
 def compact(self,n):
  order(n);self._check();a=self.budget.compact(n,self.pressure.pole)
  if set(a)!=set(('KT','KA','I')):raise ValueError('Complete original compact operator tuple required')
  return {k:EP(v) for k,v in a.items()}
 def quotient_tail(self,n):
  order(n);self._check()
  # Calling the prior rows enforces every old-field/new-pressure domain guard.
  self.current.rows(n)
  N=n+3;chart=self.state.primary.chart;lo=chart.Rmin
  W=EP(self.state.primary.covariance_R(N)['W'])
  reciprocal=EP(sum(lo**(-j-1) for j in range(N+1)))
  return prefast_quotient_rows(n,W,reciprocal,self.pressure.KT0,
   self.compact(n)['KA'],self.state.compensator.moment(self.pressure.pole),
   self.state.compensator.rho(n+1))
 def source(self,n):
  order(n);old=self.current.rows(n);tail=self.quotient_tail(n)
  whole=dict(Etheta=old['regular']['Etheta'],Ez=old['regular']['Ez']+tail['A1'])
  return dict(whole=whole,regular=old['regular'],original_prefast_A1=tail,
   exact_identity=old['exact_identity'],source_order=self.rates()['whole_source'],
   previous_200001_jet_route_required_for_this_normalized_source=False,
   actual_continuum_inputs_certified=False)
 def fast_fields(self,m):
  order(m);a=self.source(m+4);b=self.source(m+5)
  out=fast_majorants(m,a['whole']['Etheta'],a['whole']['Ez'],b['whole']['Ez'],
   self.compact(m),self.compact(m+1))
  out.update(tangential_order=self.rates()['fast_tangential'],radial_order=self.rates()['fast_radial'],
   conditional_fast_coefficient_program_materialized=True,continuum_field_enclosure_certified=False)
  return out
 def recipe(self,m):
  order(m)
  return dict(output_order=m,current_source_orders=(m+4,m+5),
   direct_prefast_A1_covariance_orders=(m+7,m+8),
   original_primary_from_initial_fields_order=m+18,
   largest_original_primary_order=m+18,
   full_pressure_order=m+6,compact_operator_orders=(m,m+1,m+4,m+5),
   obsolete_conditional_raw_A1_order=m+200009,
   obsolete_route_preserved_only_as_historical_sufficient_bound=True,
   next_obligation='recompute postfast pressure/scalars and the complete new A1 time error before the sharp slow response')
 def rates(self):
  k=Q(1,100000);D=Q(2)-4*k
  return dict(regular_source=D,original_prefast_A1_direct=Q(2),
   direct_margin_over_regular=4*k,whole_source=D,
   fast_tangential=D,fast_radial=D+1,
   actual_lower_profile_support_and_field_enclosures_required=True,
   complete_postfast_residual_or_blowup_certified=False)
