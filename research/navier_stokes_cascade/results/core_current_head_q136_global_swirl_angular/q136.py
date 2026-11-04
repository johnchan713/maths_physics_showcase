"""Conditional regular-axis/final-compact/exterior F0 angular assembly."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
outer=load('q136_outer','core_current_head_q135_exterior_swirl_angular/q135.py')
axial=load('q136_axial','core_current_head_q133_final_angular_gluing/q133.py')
exact=outer.exact;order=outer.order;lg=outer.lg
IDENTITY='same-final-leading-F0-axis-compact-through-Xstar-and-exterior-cover'
def axis(K,delta,m,binary=False):
 order(m);delta=exact(delta)
 if not 0<delta<=1:raise ValueError('Positive original axis Cauchy radius required')
 if binary:
  order(K);v=K+lg(Q(m+1))+m*lg(1/delta)
 else:
  K=exact(K)
  if K<0 or m>12:raise ValueError('Nonnegative regular-axis cap and moderate order required')
  v=K*sum(delta**(-j) for j in range(m+1))
 return dict(stage='regular-leading-axis-F0',order=m,upper=v,binary=binary,
  actual_regular_axis_K_strip_and_flat_entry_certified=False)
def compact(modulation,repair_field,xmin,xmax,identity):
 xmin,xmax=map(exact,(xmin,xmax))
 if identity!=IDENTITY or not 0<xmin<=xmax:raise ValueError('SAME final compact cover bounded away from zero required')
 axial.compact(modulation,repair_field,axial.IDENTITY) # validates FINAL stages
 P=modulation['prefactors']['E']+axial.BUMP0*repair_field['prefactor']
 return dict(stage='final-compact-F0-raw-Gevrey3',prefactor=max(Q(1),1/xmin)*P,
  base=max(modulation['base'],repair_field['base']),xmin=xmin,xmax=xmax,
  actual_whole_compact_caps_full_target_and_acceptance_certified=False)
def compact_norm(recipe,m,binary=False):
 order(m)
 if recipe.get('stage')!='final-compact-F0-raw-Gevrey3':raise ValueError('FINAL compact F0 raw recipe required')
 P,B=recipe['prefactor'],recipe['base']
 if binary:v=lg(max(Q(1),P))+lg(Q(m+1))+m*lg(max(Q(1),B))+(2*m*lg(Q(m)) if m else 0)
 else:
  if m>12:raise ValueError('Use binary budget at high order')
  v=P*sum(B**j*factorial(j)**2 for j in range(m+1))
 return dict(stage='final-compact-F0-factorial-norm',order=m,upper=v,binary=binary,xmin=recipe['xmin'],xmax=recipe['xmax'])
def exterior_norm(S0,lam,Xstar,m,binary=False):
 r=outer.exterior(S0,lam,Xstar,outer.IDENTITY)
 v=outer.norm(r,m,binary);v['Xstar']=exact(Xstar)
 return v
def join(axis_norm,compact_norm,exterior_norm,axis_end,Xstar,identity):
 axis_end,Xstar=map(exact,(axis_end,Xstar));items=(axis_norm,compact_norm,exterior_norm)
 stages=('regular-leading-axis-F0','final-compact-F0-factorial-norm','whole-exterior-leading-F0-factorial-angular-norm')
 if identity!=IDENTITY or any(v.get('stage')!=s for v,s in zip(items,stages)):raise ValueError('Same actual axis, FINAL compact and whole exterior F0 inputs required')
 m=axis_norm['order'];binary=axis_norm['binary']
 if any(v['order']!=m or v['binary'] is not binary for v in items) or Xstar<100 or axis_end<compact_norm['xmin'] or compact_norm['xmax']<Xstar or exterior_norm.get('Xstar')!=Xstar:raise ValueError('Same order/convention and overlapping full spatial cover required')
 return dict(stage='whole-leading-F0-pointwise-X-angular-norm',order=m,upper=max(v['upper'] for v in items),binary=binary,
  stronger_sum_of_global_suprema_certified=False,actual_profile_caps_cover_gates_and_scales_certified=False)
def first_inputs(F,U,m):
 order(m)
 if F.get('stage')!='whole-leading-F0-pointwise-X-angular-norm' or F.get('order')!=m or U.get('stage')!='whole-leading-pointwise-X-angular-algebra-norm' or U.get('order')!=m+2 or F.get('binary') is not U.get('binary'):raise ValueError('Q125 requires SAME-convention final F0 order m and U0 order m+2')
 return dict(F0_upper=F['upper'],U0_upper=U['upper'],F_order=m,U_order=m+2,binary=F['binary'],
  actual_profile_binding_core_collar_mixed_radial_support_or_scale_certified=False)
def record():
 return dict(status='regular-axis-final-compact-and-exterior-F0-angular-assembly',
  axis='regular F0 cap K(ell_axis) on same strip: K*sum delta^-j; never E/sqrt(2X) at zero',
  compact='PF=max(1,xmin^-1)*(Q131_E_pref+B0*Q132_field_pref), B=max(two bases); PF*sum B^j*(j!)^2',
  geometry='axis_end>=compact xmin>0, compact xmax>=Xstar>=100; all same actual profile',
  whole='intended pointwise-X factorial norm bounded by regional max; no stronger global-supremum sum claim',
  first='Q125 F0 through m and finalized U0 through m+2; no source/collar or mixed-radial proof follows',
  pressure='Q132 full baseline and all five incoming moment defects remain premises; no pressure reset',
  actual_profile_cone_pulse_covariance_residual_or_blowup_certified=False,historical_frozen_records_modified=False)
