"""Initial exact cancellation -> current Qcirc source, with unresolved A1 explicit."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
full=load('initial_cancel_full','core_full_compact_pressure_response/bounds.py')
fast=load('initial_cancel_difference','core_fast_row_gain_join/budgets.py')
EP,order,exact=full.EP,full.order,full.exact

def prefast_zero_haar_rows(n,W,graph,reciprocal,KT0):
 """ORIGINAL primary, zero means, Qcirc channel ONLY.

 Sigma and rho P disappear by exact auxiliary projection. Return the
 unprojected positive coefficient polynomial; Qcirc's factor 2 is kept
 separately and is already present in the original A1 raw recipe.
 """
 order(n);W,g,r,kp=map(EP,(W,graph,reciprocal,KT0));d=n+4
 gr=(d*g+d+2*r)*W;pm=kp*gr
 return dict(Etheta=(d*g+d+2*r)*W,Ez=(d*g+r)*W+d*(W+pm),
  Qcirc_norm_factor=2,covariance_input_order=n+3)

class CurrentFastSource:
 """Original initial fields are read from their actual lower norm recipes.

 No fine current-row callback is accepted. The unresolved original A1
 tail is kept as a separate exact raw-kernel/source obligation, never
 forwarded as if it vanished or already had the proposed fine rate.
 """
 def __init__(self,pressure,axial_gap):
  if not isinstance(pressure,full.FullCompactPressure):raise TypeError('Full auxiliary compact-pressure program required')
  self.pressure=pressure;self.state=pressure.state;self.gap1=exact(axial_gap)
  if not 0<self.gap1<=self.state.target.Llo/2:raise ValueError('Actual ORIGINAL axial A1 cutoff gap required')
  h=self.state.haar
  if not isinstance(h.initial,full.recompute.projector.haar.sharp.SharpInitial):raise TypeError('Original sharp initial field recipe required')
  self.initial=h.initial;self.budget=h.initial.budget
  if self.budget.primary is not self.state.primary or self.budget.pole!=h.pole or self.state.target.pole!=2*h.pole:raise ValueError('Initial fields, Haar rows and signed update must share the same prefix')
  self.initial_pole=self.budget.pole
 def tail_recipe(self,n):
  order(n);r=fast.flat_recipe(n);J=r['raw_source_order']
  return dict(r,stage='ORIGINAL initial A1 Qcirc Ez_prefast',axial_cutoff_gap=self.gap1,
   original_width=self.state.target.width,
   original_factors=['original U integration length','N_n(k1 chi_U)','N_J(f1)',
    'ORIGINAL primary prefast Ez bound at J','Eedge(source pole,original width)',
    'inverse weight floor on the ORIGINAL axial cutoff support'],
   floor_factor='(1+initial_compact_I(0)*S)^integration_order',
   source_program=dict(method='prefast_zero_haar_rows on ORIGINAL complete primary and FULL shifted T0',
    row_order=J,covariance_order=J+3,full_pressure_order=J+1,independent_fine_source_callback=False),
   minimum_primary_source_pole=12*(J+3)+18,
   existing_initial_pole=self.budget.pole,
   existing_pole_covers_this_source=self.budget.pole>=12*(J+3)+18,
   time_floor='(1+initial_compact_I(0)*S)^integration_order; no new temporal inverse',
   raw_cutoff_and_multiplier_bounds_still_required=True,actual_tail_polynomial_materialized=False)
 def rows(self,n):
  order(n);s=self.state;s.check();N=n+3
  if s.haar.initial is not self.initial or self.initial.budget is not self.budget or self.budget.pole!=self.initial_pole or s.haar.pole!=self.initial_pole:raise ValueError('Original initial prefix or proof parameter changed')
  if n>s.haar.max_output or self.pressure.order<n+1:raise ValueError('Current source/pressure derivative ledger too short')
  if n>self.initial.max_output or self.budget.pole<12*(n+13)+18:raise ValueError('Original initial field source/covariance ledger too short')
  if self.pressure.pole<2*self.budget.pole:raise ValueError('Initial fast quadratic/pressure domain not covered')
  f={k:EP(v) for k,v in self.budget._fast(N).items()}
  if set(f)!={'x','u','z'}:raise ValueError('Actual initial fast radial/angular/axial fields required')
  v=self.initial.slow_fields(N,True);slow=dict(x=EP(v['X']),u=EP(v['U']),z=EP(v['Z']))
  base=self.budget.fixed(N,self.pressure.pole)
  if not {'b','V','G','c'}<=set(base):raise ValueError('Original base and frozen slow-time coefficients required')
  chart=s.primary.chart;lo=chart.Rmin
  fixed={k:EP(base[k]) for k in ('b','V','G','c')}
  fixed.update(graph=EP(2+6*chart.B*sum((j+1)*lo**(-j) for j in range(N+1))),
   r=EP(sum(lo**(-j-1) for j in range(N+1))),rho=s.compensator.rho(n+1),
   Kp=self.pressure.KT0,Kmom=EP(s.compensator.moment(self.pressure.pole)))
  zero={k:EP(0) for k in ('x','u','z')}
  df=fast.difference(n,zero,f,fixed)
  ds=fast.difference(n,f,slow,fixed)
  dw=self.pressure.covariance_response(n)
  # Qcirc kills the original prefast Haar row. A1 maps zero Haar to zero
  # Haar exactly, so its separate tail needs no additional Qcirc factor.
  regular={k:2*(df[k]+ds[k]+dw[k]) for k in ('Etheta','Ez')}
  return dict(regular=regular,original_fast_difference=df,original_slow_difference=ds,
   complete_signed_covariance_difference=dw,axial_A1_tail=self.tail_recipe(n),
   exact_identity='Qcirc EpostSigned_regular=(0,A1 Qcirc Ez_prefast)+Qcirc(Rinitial_fast+Rinitial_slow+Rcov)',
   regular_angular_source_complete_under_original_input_premises=True,
   axial_source_regular_plus_unmaterialized_A1=True,actual_whole_fine_source_certified=False)
 def fast_recipe(self,m):
  order(m)
  return dict(stage='new fast after FIRST real signed update',output_order=m,
   tangential_source_order=m+4,radial_source_order=m+5,
   maximum_new_signed_covariance_order=m+8,maximum_new_signed_target_order=m+9,
   original_primary_from_initial_fields_order=m+18,
   original_primary_from_A1_tail_order=self.tail_recipe(m+5)['source_program']['covariance_order'],
   minimum_recomputed_pressure_max_output=m+6,sufficient_profile_interface_order=m+10,
   angular_source='CurrentFastSource.rows(m+4).regular.Etheta; original A1 tail is axial only',
   axial_source='CurrentFastSource.rows(...).regular.Ez PLUS original A1 tail',
   actual_new_axial_velocity='gd-A1 gd; never gd alone',
   actual_fresh_fast_fields_materialized=False)
 def inverse_source(self,*args,**kw):
  raise NotImplementedError('Original A1 raw-kernel/profile/source polynomial and whole-domain bounds must be instantiated before a materialized fine axial inverse input')
 def rates(self):
  k=Q(1,100000);B=1-k;alpha=2-3*k;D=alpha-k
  return dict(initial_fast_difference=alpha,initial_slow_difference=2-2*k,changed_covariance=D,
   regular_current_source=D,original_A1_coarse=B,original_A1_conditional_sharp=B+99999*k,
   whole_source_without_raw_A1_gain=B,whole_source_if_original_A1_obligations_hold=D,
   actual_fine_source_or_mean_closure_certified=False)
