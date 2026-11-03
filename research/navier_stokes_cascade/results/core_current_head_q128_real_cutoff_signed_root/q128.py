"""Real Gevrey-3 angular cutoff/root bounds; full pre-modulation memory is an input."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
current=load('q128_current','core_current_head_q127_annular_radius_axis/q127.py')
signed=load('q128_signed','core_current_head_q65_signed_root_radius/q65.py')
order=current.order;exact=current.exact;lg=current.lg
FIELDS={'E','U','EX','UX','invE'};TRACES={'M','I','J','S','Pi'}
IDENTITY='same-real-variance-cutoff-and-signed-root-with-full-forward-memory'

def state_caps(h,x0,x1,rho,fields,traces):
 h,x0,x1,rho=map(exact,(h,x0,x1,rho))
 if not 0<h<Q(1,100) or not 0<x0<=x1 or not 0<rho<=Q(1,2):raise ValueError('Original positive compact interval and complex eta radius required')
 if set(fields)!=FIELDS or set(traces)!=TRACES:raise ValueError('All field/inverse and all five incoming complex trace caps required')
 f={k:exact(v) for k,v in fields.items()};t={k:exact(v) for k,v in traces.items()}
 if min(f.values())<0 or min(t.values())<0 or min(f['E'],f['invE'])<=0:raise ValueError('Nonnegative caps and positive same-field E/inverse bounds required')
 E,U,EX,UX,inv=f['E'],f['U'],f['EX'],f['UX'],f['invE'];length=x1-x0;rcap=max(Q(1),2*x1)
 densities=dict(M=U,I=rcap*E,J=rcap*E*U,S=U*U+E*E/2,Pi=E*E/(2*x0))
 mem={k:t[k]+length*densities[k] for k in TRACES}
 D=Q(1,2)-h;A=Q(1,2)+h;M,S,Pi=mem['M'],mem['S'],mem['Pi']
 W=x1+(4*D+10/rho)*M
 p=(W*U+D*(1+4/rho)*M+(8*h+10/rho)*S+x1*(8*A+10/rho)*Pi)*inv/(1-8*h)
 return dict(stage='PREMOD-state',rho=rho/2,a=1+2*x1*EX*inv,b=2*x1*UX*inv,p=p,
             primitives=mem,densities=densities,full_incoming_pressure_required=True,
             actual_uniform_complex_fields_and_traces_certified=False)

def parameters(state,afloor,delta,G,E_exp,d0,identity):
 if identity!=IDENTITY or state.get('stage')!='PREMOD-state':raise ValueError('Same real cutoff/root applied to the complete PREMOD state required')
 rho,A,B,P=map(exact,(state['rho'],state['a'],state['b'],state['p']))
 afloor,delta,G,E_exp,d0=map(exact,(afloor,delta,G,E_exp,d0))
 if not 0<rho<=Q(1,2) or not 0<afloor<=A or min(B,P)<0 or delta<=0:raise ValueError('Original real a floor, positive variance gap and complex state caps required')
 ra=min(rho/2,rho*afloor/(4*A));V=A+2*B*B/afloor;gap=delta/4;Hstar=2+delta/2
 rh=min(ra/2,ra*gap/(4*V));C=2*max(Q(1),Hstar/afloor)
 Lz=(1+64*V/delta)/ra;Ls=Lz+1/rh;Ms=729*C
 inv=signed.radii(G,E_exp,d0)
 Lmu=Ls*(1+2*(P+Ms)/inv['sigma'])
 return dict(stage='real-Gevrey3-cutoff-root',a_radius=ra,v_cap=V,radicand_radius=rh,
             amplitude_prefactor=Ms,amplitude_base=Ls,mu_prefactor=inv['a'],mu_base=Lmu,
             scalar_inverse_radius=inv['sigma'],real_mu_value_cap=G,factorial_power=3,
             slow_analytic_extension_claimed=False,
             actual_real_floors_value_caps_or_E_exp_bound_certified=False,
             actual_final_compact_profile_or_scale_certified=False)

def derivative_bound(recipe,m,kind,binary=False):
 order(m)
 if recipe.get('stage')!='real-Gevrey3-cutoff-root' or kind not in ('amplitude','mu'):raise ValueError('Cutoff amplitude or same signed-root recipe required')
 if kind=='mu' and m==0:
  pref=recipe['real_mu_value_cap'];base=Q(1)
 else:pref=recipe[kind+'_prefactor'];base=recipe[kind+'_base']
 if binary:
  value=lg(max(Q(1),pref))+m*lg(max(Q(1),base))+(3*m*lg(Q(m)) if m else 0)
 else:
  if m>12:raise ValueError('Use the binary factorial budget at high order')
  value=pref*base**m*factorial(m)**3
 return dict(order=m,kind=kind,upper=value,binary=binary,factorial_power=3,
             slow_analytic_extension_claimed=False,actual_input_bindings_certified=False)

def record():
 return dict(status='real-cutoff-and-signed-root-Gevrey3-angular-bound',
  amplitude='|s^(m)| <= Ms*Ls^m*(m!)^3',
  signed_positive_jets='|mu^(m)| <= a_Psi*[Ls*(1+2*(P+Ms)/sigma_Psi)]^m*(m!)^3',
  signed_value='original G bound retained at m=0',
  full_forward_pressure_terms='X*(4*Aphys*eta*Pi-(1-eta^2)*Pi_eta)/(L*E)',
  scalar_inverse_E_exp_outward_bound_is_a_premise=True,
  slow_analytic_extension_claimed=False,
  actual_final_compact_profile_scale_cone_pulse_residual_or_blowup_certified=False)
