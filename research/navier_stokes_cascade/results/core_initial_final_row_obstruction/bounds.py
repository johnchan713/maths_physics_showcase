"""Exact first-real final-row join and its sharp zero-margin obstruction."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

previous=load('initial_final_previous','core_initial_postfast_refresh/bounds.py')
closure=load('initial_final_closure','core_final_mean_row_closure/budgets.py')
relaxation=load('initial_final_relaxation','post_signed_mean_relaxation/relax.py')
EP,order,exact=previous.EP,previous.order,previous.exact
FIELDS=('x','u','z')

def a1_phi_recipe(n,*,geometry,pole,kappa=Q(1,100000)):
 """Raw compact recipe for A1 Dz Q Phi at its SAME post-signed state.

 The generic raw source index belongs to Dz Phi.  Phi itself needs one
 additional derivative, which is why its index is n+2p+4.
 """
 order(n);q=previous.a1_time_error_request(n,kappa);p=q['integration_order']
 raw=closure.oscillatory_recipe(n,p,e=1,source='Dz Qcirc Phi_postsigned',
  source_rate=q['source_order'],target_rate=q['target_order'],kappa=kappa,
  geometry=geometry,pole=pole)
 raw.update(raw_Dz_Phi_order=raw['raw_source_order'],raw_Phi_order=raw['raw_source_order']+1,
  identity=q['identity'],Phi=q['Phi'],source_stage=q['source_stage'],
  strict_margin=q['strict_margin'],
  postfast_or_postslow_pressure_substituted=False,
  potential_dependencies={k:q[k] for k in ('current_fields_for_potential','current_Ez_for_potential',
   'coarse_pressure_for_potential','conservative_pressure_field_recipe',
   'current_Etheta_for_fine_pressure','fixed_G_for_GR')},
  actual_Phi_profile_support_and_high_jet_bounds_materialized=False)
 return raw

def final_row_envelope(n,post_signed_Haar,A1_tail,old_full,fast_increment,slow_increment,fixed):
 """Exact positive envelope for Haar(EpostSigned)+(0,A1Fz)+Rf+Rs."""
 order(n)
 if set(post_signed_Haar)!=set(('Etheta','Ez')):raise ValueError('Both complete post-signed Haar rows required')
 if any(set(v)!=set(FIELDS) for v in (old_full,fast_increment,slow_increment)):
  raise ValueError('Three complete old, fast and slow mean components required')
 f=previous.previous.prior.fast.difference(n,old_full,fast_increment,fixed)
 after={k:EP(old_full[k])+EP(fast_increment[k]) for k in FIELDS}
 s=previous.previous.prior.fast.difference(n,after,slow_increment,fixed)
 return dict(Etheta=EP(post_signed_Haar['Etheta'])+f['Etheta']+s['Etheta'],
  Ez=EP(post_signed_Haar['Ez'])+EP(A1_tail)+f['Ez']+s['Ez'],
  post_signed_Haar={k:EP(v) for k,v in post_signed_Haar.items()},A1=EP(A1_tail),
  fast_difference=f,slow_difference=s,
  exact_identity='Efinal=Haar(Epostsigned)+(0,A1 Qcirc Ez_postsigned)+Rf+Rs',
  extra_Qcirc_factor_on_full_differences=False)

def complete_final_pressure(postfast_pressure,slow_pressure_difference):
 """Append the sharp slow pressure difference to all postfast memory."""
 if not isinstance(postfast_pressure,dict) or 'pm' not in postfast_pressure or 'parts' not in postfast_pressure:
  raise ValueError('Complete postfast pressure decomposition required')
 parts=dict(postfast_pressure['parts']);parts['sharp_slow']=EP(slow_pressure_difference)
 return dict(parts=parts,pm=EP(postfast_pressure['pm'])+parts['sharp_slow'],
  complete_pressure_memory_retained=True,additive_pressure_gauge_separate=True)

def rate_ledger(kappa=Q(1,100000)):
 k=exact(kappa)
 if not 0<k<=Q(1,100):raise ValueError('Retained positive small-kappa regime required')
 B=1-k;alpha=2-3*k;D=alpha-k;nu=3-6*k
 source=dict(initial_fast_difference=alpha,initial_slow_difference=2-2*k,
  complete_signed_covariance_response=D,corrected_prefast_A1=Q(2))
 fast=dict(slow_time=D+1,base_radial_flux=D+1-k,radial_viscosity=D+1-2*k,
  old_mean_transport=B+D+1-k,changed_axial_pressure=D+1,new_mean_square=2*D+1-k)
 slow=dict(slow_time=alpha+1,base_radial_flux=alpha+1,radial_viscosity=alpha+1,
  old_mean_transport=B+alpha+1-k,changed_axial_pressure=alpha+1,new_mean_square=2*alpha+1)
 return dict(B=B,alpha=alpha,D=D,nu=nu,source=source,fast=fast,slow=slow,
  source_margins={name:value-D for name,value in source.items()},
  fast_margins={name:value-nu for name,value in fast.items()},
  slow_margins={name:value-nu for name,value in slow.items()},
  whole_source=min(source.values()),fast_difference=min(fast.values()),slow_difference=min(slow.values()),
  unique_zero_margin_fast_groups=tuple(name for name,value in fast.items() if value==nu),
  source_zero_margin_groups=tuple(name for name,value in source.items() if value==D))

def strict_fast_closure(delta_source=0,delta_viscosity=0,kappa=Q(1,100000)):
 """Necessary inequality inside the inherited positive-majorant route."""
 a,b=map(exact,(delta_source,delta_viscosity));r=rate_ledger(kappa)
 if min(a,b)<0:raise ValueError('Only proved nonnegative gains may be inserted')
 achieved=r['fast']['radial_viscosity']+a+b
 return dict(delta_source=a,delta_viscosity=b,achieved=achieved,target=r['nu'],
  strict_margin=achieved-r['nu'],strictly_closes=achieved>r['nu'],
  required_inequality='delta_source+delta_viscosity>0',
  alternatives=('improve the complete signed-covariance source above D',
   'prove an exact cancellation or positive epsilon gain in the radial viscosity contribution'),
  no_positive_gain_present_in_inherited_bounds=a+b==0)

class InitialFinalRowObstruction:
 def __init__(self,postfast):
  if not isinstance(postfast,previous.InitialPostfastRefresh):
   raise TypeError('Complete initial postfast refresh required')
  self.postfast=postfast
 def final_rows(self,n,post_signed_Haar,A1_tail):
  order(n);N=n+3
  return final_row_envelope(n,post_signed_Haar,A1_tail,self.postfast.old_full(N),
   self.postfast.new_fast(N),self.postfast.slow_fields(N),self.postfast._fixed(n))
 def final_pressure(self,n,post_signed_Haar,A1_tail):
  rows=self.final_rows(n,post_signed_Haar,A1_tail)
  return complete_final_pressure(self.postfast.pressure_envelope(n),rows['slow_difference']['pm'])
 def recipe(self,n,*,geometry,pole):
  return dict(output_order=n,final_identity=final_row_envelope.__doc__,
   A1=a1_phi_recipe(n,geometry=geometry,pole=pole),rates=rate_ledger(),
   strict_obstruction=strict_fast_closure(),
   next_obligation='prove a positive source/viscosity gain or exact cancellation, while validating the deferred same-state Phi high jets and support')
