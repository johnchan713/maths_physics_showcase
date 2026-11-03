"""Earlier annular PREMOD radius and the same leading axis angular cover."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
current=load('q127_current','core_current_head_q126_exterior_axial_radius/q126.py')
old=load('q127_earlier_root','core_current_head_q73_annular_target_all_orders/q73.py')
order=current.order;exact=current.exact;lg=current.current.lg
PREMOD_ID='same-earlier-eta-dependent-annular-PREMOD-root'
AXIS_ID='same-leading-core-through-flat-modulation-entry'
FINAL_STAGE='final-leading-compact-modulation-and-repair'
BETA=Q(2000);QCAP=Q(200);D=Q(6,10**13);A=1/(16*BETA*QCAP)

def radius(H,Z,identity):
 H,Z=map(exact,(H,Z))
 if identity!=PREMOD_ID or not 0<H<=1 or Z<0:raise ValueError('Same earlier target strip, nonnegative target cap and PREMOD identity required')
 S=2*Z/H+QCAP*D*D
 sigma=min(H/2,Q(1,2),A/(8*BETA*S))
 linear=2*BETA*QCAP*sigma*D
 checks=dict(same_root_inside=2*D<A,constant_term=BETA*S*sigma<=A/8,
             image=BETA*S*sigma+linear*A+BETA*QCAP*A*A<A/2,
             lipschitz=linear+2*BETA*QCAP*A<Q(1,4),premod_field_cap=8+128*(D+A)<16,relative_swirl_nonzero=192*(D+A)<Q(1,2))
 assert all(checks.values())
 return dict(sigma=sigma,increment_radius=A,real_value_cap=D,checks=checks,
             actual_target_and_same_root_contracts_certified=False)

def inverse_radius_bits(H_inverse_bits,Z_bits):
 order(H_inverse_bits);order(Z_bits)
 # 8*beta/A <2^37, S<=2^(2+Z_bits+H_inverse_bits).
 return 39+H_inverse_bits+Z_bits

def premod_U(m,ua,H,Z,identity,binary=False):
 order(m)
 if binary:
  if identity!=PREMOD_ID:raise ValueError('Earlier PREMOD root required')
  order(ua);inv=inverse_radius_bits(H,Z)
  value=1+max(ua,4)+lg(Q(m+1))+m*inv
 else:
  ua=exact(ua)
  if ua<0 or m>40:raise ValueError('Nonnegative moderate PREMOD value cap; use bits at high order')
  r=radius(H,Z,identity);value=(ua+16)*sum(r['sigma']**(-j) for j in range(m+1))
 return dict(stage='PREMOD',order=m,upper=value,binary=binary,
             actual_prefix_target_or_final_modulated_norm_certified=False)

def axis_U(m,axis_value,delta,identity,binary=False):
 order(m);delta=exact(delta)
 if identity!=AXIS_ID or not 0<delta<=1:raise ValueError('Same leading axis cover and positive original Cauchy radius required')
 if binary:
  order(axis_value);value=axis_value+lg(Q(m+1))+m*lg(1/delta)
 else:
  axis_value=exact(axis_value)
  if axis_value<0 or m>40:raise ValueError('Nonnegative moderate axis bound; use bits at high order')
  value=axis_value*sum(delta**(-j) for j in range(m+1))
 return dict(stage='leading-axis',order=m,upper=value,binary=binary,
             actual_original_axis_cap_and_cover_certified=False)

def leading_U(m,axis_value,delta,final_compact,axis_identity,pulse_identity,binary=False):
 a=axis_U(m,axis_value,delta,axis_identity,binary)
 if not isinstance(final_compact,dict) or final_compact.get('stage')!=FINAL_STAGE or final_compact.get('order')!=m or final_compact.get('binary') is not binary:
  raise ValueError('A separately supplied FINAL compact norm at this same order/convention is required; PREMOD is not that input')
 g=current.global_U(m,a['upper'],final_compact['upper'],pulse_identity,binary)
 return dict(axis=a,global_U=g,actual_final_compact_norm_and_scale_certified=False)

def record():
 return dict(status='same-earlier-annular-radius-and-axis-angular-cover',
  increment_radius='a=1/(16*2000*200)',real_accepted_value_cap='6e-13',
  radius='min(H/2,1/2,a/[8*2000*(2Z/H+200*d^2)])',
  coefficient_positive_jet='a*sigma^-k for k>=1; original real value bound d retained',
  earlier_annular_swirl='|E|<=4P, |1/E|<=20/P on -6<=y<=0',
  PREMOD_U='(ua+16)*sum_(j=0)^m sigma^-j',axis_U='K_axis_value*sum_(j=0)^m delta^-j',
  binary_inverse_radius='39+H_inverse_bits+Z_bits',
  actual_final_compact_profile_scale_cone_pulse_residual_or_blowup_certified=False)
