"""Fixed angular radius for the SAME exterior axial root, under its original contracts."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
current=load('q126_current','core_current_head_q125_uniform_angular_first_seed/q125.py')
pulse=load('q126_pulse','core_current_head_q77_axial_profile_pulse/q77.py')
order=current.order;exact=current.exact
IDENTITY='same-exterior-positive-energy-root-and-preheat-future-mass'

def value_checks():
 L=Q(1,2**256);tL=Q(1,2**250);r=2**19*L**30;Z=2**15*L**30
 return dict(exponential_upper=Q(65,24)+Q(1,100)<Q(11,4),
  target_disk_constant=10*Q(11,4)**8<2**15,
  angular_ball=r<=Q(1,1024),same_accepted_root=7000*L**30<=r,
  angular_image=5*Z+160*r*r<r,angular_lipschitz=320*r<Q(1,2),
  future_mass_factor=(1+128*r)**2<2,
  affine_end_energy=1960*L**40<64,
  energy_coefficient=Q(123,500)+40*L**41<Q(1,4),
  energy_forcing=2**15*tL<Q(1,1000),
  energy_slope_change=1024*L**41<Q(1,1000),
  energy_map_image=3*(Q(1,1000)+Q(1,1000)*Q(1,10)+Q(1,4)*Q(1,100))<Q(1,10),
  energy_map_lipschitz=3*(Q(1,1000)+Q(1,20))<Q(1,2),
  energy_constant_negative=2**14*tL<(1-Q(1,2**26))/4,
  pulse_complex_cap=39+1792*L**20<40)

def axial_angular(m,identity,binary=False):
 order(m)
 if identity!=IDENTITY:raise ValueError('The same exterior pulse/pre-heat root contract is required')
 if binary:value=m+7
 else:
  if m>40:raise ValueError('Use the binary bound at high angular order')
  value=40*(2**(m+1)-1)
 return dict(order=m,upper=value,binary=binary,angular_disk_radius=Q(1,2),
             actual_parameter_root_and_stage_contracts_certified=False,
             heat_tail_complex_extension_claimed=False,physical_wave_bound_certified=False)

def global_U(m,axis,compact,identity,binary=False):
 a=axial_angular(m,identity,binary)
 if binary:
  order(axis);order(compact);value=max(axis,compact,3,a['upper'])
 else:
  axis,compact=map(exact,(axis,compact))
  if min(axis,compact)<0:raise ValueError('Original axis and compact-region norms required')
  value=max(axis,compact,Q(8),Q(a['upper']))
 return dict(order=m,upper=value,binary=binary,exterior=a,
             actual_axis_compact_matching_or_support_inputs_certified=False)

def current_first(h,m,nu,Xcut,Xvel,Xv,axis_U,compact_U,F0,delta,core,collar,
                  pulse_identity,seed_identity,binary=False):
 U=global_U(m+2,axis_U,compact_U,pulse_identity,binary)
 first=current.first_discrepancies(h,m,nu,Xcut,Xvel,Xv,U['upper'],F0,delta,core,collar,seed_identity,binary)
 return dict(global_U=U,first=first,actual_leading_profile_or_scale_certified=False)

def record():
 checks=value_checks()
 assert all(checks.values())
 return dict(status='fixed-exterior-axial-angular-radius-and-current-first-join',
  angular_radius='1/2 at every real eta center',angular_root_ball='2^19*lambda^30',
  preheat_future_mass_complex_cap='2^14*T',strength_complex_cap='3/2',
  exterior_U_complex_cap=40,exterior_U_angular_norm='40*(2^(m+1)-1)',
  exterior_U_binary_upper='m+7',value_checks=checks,
  actual_axis_compact_profile_scale_cone_pulse_residual_or_blowup_certified=False,
  heat_tail_complex_extension_claimed=False)
