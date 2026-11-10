"""Value-based mixed Gevrey-3 recurrence on Q71's post-taper continuation."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
previous=load('q143_previous','core_current_head_q142_short_taper_order_zero/q142.py')
mixed=previous.previous
exact=mixed.exact
IDENTITY='same-Q142-endpoint-and-Q71-post-taper-prescription'

def endpoint_signature(short):
 if short.get('stage')!='Q141-short-taper-with-certified-order-zero-inputs':
  raise ValueError('Same certified Q142 endpoint required')
 r=short['recipe']
 rows=[]
 for group,names in ((r['reference'],('F','invF','U')),(r['sources']['moments'],('M','I','J','S','Pi')),(r['state']['moments'],('M','I','J','S','Pi'))):
  rows.extend(tuple(group[n][k] for k in ('P','A','B')) for n in names)
 return tuple(rows)+(short['a_floor'],short['a_cap'])

def reference_sources(short,h,t,Xi,Lambda,length,Benv,identity):
 """Reintegrate Q141 reference memory, then use Q71's exact modes."""
 h,t,Xi,Lambda,length,Benv=map(exact,(h,t,Xi,Lambda,length,Benv))
 if identity!=IDENTITY or not 0<h<Q(1,2) or min(t,Xi,Lambda,length,Benv)<=0 or Xi>110 or length>Benv or length.denominator!=1 or Xi*Lambda/4>2**length.numerator:
  raise ValueError('Same positive Q142/Q71 continuation data required')
 if short.get('a_floor')!=2*t/(5*Lambda) or short.get('a_cap')!=Benv/100:
  raise ValueError('Same Q142 endpoint floor, Lambda and Benv required')
 r=short['recipe'];ref=r['reference'];memory=r['sources']['moments']
 Be=max(ref[n]['B'] for n in ('F','invF','U'))
 flat=dict(stage='Q72-value-based-field-recurrence',
  F=mixed.env(ref['F']['P'],1,Be,'constant-in-y-Gevrey3'),
  invF=mixed.env(ref['invF']['P'],1,Be,'constant-in-y-Gevrey3'),
  U=mixed.env(ref['U']['P'],1,Be,'constant-in-y-Gevrey3'))
 T=max(v['P'] for v in memory.values());Bt=max(v['B'] for v in memory.values())
 traces=mixed.trace_rows(T,Bt);geo=mixed.geometry(Xi,Lambda/4,15,Lambda)
 raw=mixed.sources(flat,traces,h,geo,length,mixed.IDENTITY)
 # Q71's exact three-point inversion has row sums <=10 for p1r and
 # <=19 for X*nsr; their positive y derivatives cost at most 1^i,2^i.
 return dict(stage='Q71-endpoint-derived-three-mode-Gevrey3',
  p1=mixed.env(10*raw['p1']['P'],1,raw['p1']['B']),
  Xns=mixed.env(19*raw['Xns']['P'],2,raw['Xns']['B']),raw=raw,
  endpoint_signature=endpoint_signature(short),parameters=(h,t,Xi,Lambda,length,Benv),
  all_five_reference_rows_reintegrated=True,
  full_reference_I_J_Pi_and_Pi_eta_memory_retained=True)

def post_taper(short,source,h,t,Xi,Lambda,length,Benv,identity):
 if identity!=IDENTITY or short.get('stage')!='Q141-short-taper-with-certified-order-zero-inputs' or source.get('stage')!='Q71-endpoint-derived-three-mode-Gevrey3':
  raise ValueError('Same certified short-taper endpoint and Q71 source modes required')
 h,t,Xi,Lambda,length,Benv=map(exact,(h,t,Xi,Lambda,length,Benv))
 if not 0<h<Q(1,2) or min(t,Xi,Lambda,length,Benv)<=0 or Xi>110 or length>Benv or length.denominator!=1 or Xi*Lambda/4>2**length.numerator or t*Benv*Benv>Q(1,10):
  raise ValueError('Selected post-taper interval and t*B^2 smallness required')
 # Q142 did not copy its input tuple into the returned record, but it did
 # freeze these two identities.  They recover the same Lambda and Benv and
 # prevent joining a Q141 endpoint to a different post-taper geometry.
 if short.get('a_floor')!=2*t/(5*Lambda) or short.get('a_cap')!=Benv/100:
  raise ValueError('Same Q142 endpoint floor, Lambda and Benv required')
 if source.get('endpoint_signature')!=endpoint_signature(short) or source.get('parameters')!=(h,t,Xi,Lambda,length,Benv):
  raise ValueError('Reference sources and actual continuation must share one endpoint and parameter tuple')
 r=short['recipe'];state0=r['state'];actual0=r['actual']
 step=mixed.switch(t,mixed.IDENTITY)
 # beta, theta and 1-theta have the same fixed Gevrey-3 step envelope.
 a_src=mixed.add(mixed.scale(mixed.mul(step,source['p1']),t),mixed.scale(step,Q(4,5)))
 u_src=mixed.scale(mixed.mul(step,source['Xns']),t/2)
 Btrace=max(actual0['F']['B'],actual0['invF']['B'],actual0['U']['B'])
 Ftrace=mixed.env(actual0['F']['P'],1,Btrace,'constant-in-y-Gevrey3')
 invFtrace=mixed.env(actual0['invF']['P'],1,Btrace,'constant-in-y-Gevrey3')
 Utrace=mixed.env(actual0['U']['P'],1,Btrace,'constant-in-y-Gevrey3')
 zero=mixed.env(0,1,max(a_src['B'],u_src['B']),'constant-in-y-Gevrey3')
 L=mixed.integrate(mixed.scale(a_src,Q(1,2)),zero,length)
 # F/F(B)=exp(-int a/2).  Before X=100 the exponent is at most
 # t B^2/200; on [100,110], a<=4/5 and log(11/10)<1/10.
 # The reciprocal multiplier is therefore <exp(1/20)<3.
 eL=mixed.exp_comp(L,3)
 field=dict(stage='Q72-value-based-field-recurrence',F=mixed.mul(Ftrace,eL),
  invF=mixed.mul(invFtrace,eL),U=mixed.integrate(u_src,Utrace,length),
  log_increment=L,exp_value_cap=Q(3),order_dependent_norm_exponentiated=False)
 # The ACTUAL five endpoint rows at B, not reference rows or reset values.
 T=max(v['P'] for v in state0['moments'].values())
 Bt=max(v['B'] for v in state0['moments'].values())
 traces=mixed.trace_rows(T,Bt)
 geo=mixed.geometry(Xi,Lambda/4,15,Lambda)
 out=mixed.state(field,traces,h,geo,length,mixed.IDENTITY)
 afloor=2*t/(5*Lambda)
 return dict(stage='Q71-post-taper-value-based-mixed-Gevrey3',field=field,state=out,
  reference_sources=source,incoming_actual_moments=traces,a_source=a_src,U_source=u_src,
  H=out['H'],Ay=out['Ay'],Be=out['Be'],a_floor=afloor,a_cap=Benv,
  all_five_actual_incoming_rows_retained=True,full_p2_Pi_and_Pi_eta_retained=True,
  inherited_exp_Tm_squared_cap_used=False,selected_order_zero_floor_and_exp_value_cap_certified=True,
  shape_and_annular_regional_four_field_maximum_certified=False)

def old_q73_log_lower(t,m):
 t=exact(t)
 if t<=0 or not isinstance(m,int) or isinstance(m,bool) or m<=0:raise ValueError('Positive width and order required')
 T=162*3260**m*factorial(m)**2/t**m
 # Q73 H contains (tau+1) T_m^2 Q1/2; use positive factors one.
 v=T*T/2
 return v.numerator//v.denominator

def bound(recipe,kind,i,k,binary=False):
 if recipe.get('stage')!='Q71-post-taper-value-based-mixed-Gevrey3' or kind not in ('a','b','p','E'):
  raise ValueError('Same post-taper four-field recipe required')
 return mixed.bound(recipe['state'][kind],i,k,binary)

def record():
 return dict(status='Q71-post-taper-value-based-mixed-Gevrey3-recurrence',
  correction='Q73 finite-order exp(T_m^2) does not imply fixed Gevrey constants; exponentiate only the actual order-zero field multiplier',
  reference='reintegrate the Q141 reference endpoint and all five rows on the flat continuation, then use Q71 exact three-mode radial spans for p1r and X*nsr',
  actual='same beta/theta/kappa prescriptions build F,U; all five actual moment rows enter from the Q141 endpoint and are reintegrated',
  floor='a>=2*C^-20/(5*Lambda) on B<=X<=110; final interpolation is between positive endpoints',
  output='fixed mixed Gevrey-3 H,Ay,Be and selected a floor for actual post-taper a,b,p2,E',
  blocking_region='shape and earlier annular restoration still need a common four-field mixed recurrence and gate binding before the Q140-Q138 compact join',
  actual_complete_profile_pulse_covariance_residual_or_blowup_certified=False,historical_frozen_records_modified=False)
