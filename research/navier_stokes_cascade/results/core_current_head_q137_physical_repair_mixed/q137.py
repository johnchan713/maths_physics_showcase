"""Closed ordinary-X/eta Gevrey-3 bounds for the SAME leading repair only."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
previous=load('q137_previous','core_current_head_q136_global_swirl_angular/q136.py')
tuple_provider=load('q137_tuple','core_current_head_q121_all_order_radial_tuple/q121.py')
exact=previous.exact;order=previous.order;lg=previous.lg
IDENTITY='same-leading-five-row-fixed-unit-mass-bump-repair-and-full-pressure-target'
D=Q(32600);A=Q(162)*D;G=8*(1+D)
def mixed(field,xmin,identity):
 xmin=exact(xmin)
 if identity!=IDENTITY or field.get('stage')!='same-fixed-bump-repair-field-Gevrey3' or xmin<=0:raise ValueError('Same Q132 leading repair and positive lower support radius required')
 P,K=map(exact,(field['prefactor'],field['base']))
 if P<0 or K<1:raise ValueError('Nonnegative repair prefactor and angular base >=1 required')
 return dict(stage='same-leading-physical-repair-mixed-Gevrey3',
  prefactors=dict(E=A*P,U=A*P,F=A*P*max(Q(1),1/xmin)),
  radial_bases=dict(E=G/xmin,U=G/xmin,F=(G+1)/xmin),angular_base=K,xmin=xmin,
  actual_full_target_root_support_and_positive_field_certified=False,
  complete_modulated_profile_or_pressure_support_certified=False)
def bound(recipe,kind,r,k,binary=False):
 order(r);order(k)
 if recipe.get('stage')!='same-leading-physical-repair-mixed-Gevrey3' or kind not in ('E','F','U'):raise ValueError('Same leading repair component required')
 P,R,K=recipe['prefactors'][kind],recipe['radial_bases'][kind],recipe['angular_base']
 if binary:v=lg(max(Q(1),P))+r*lg(max(Q(1),R))+k*lg(max(Q(1),K))+(3*r*lg(Q(r)) if r else 0)+(3*k*lg(Q(k)) if k else 0)
 else:
  if max(r,k)>12:raise ValueError('Use binary budget for large orders')
  v=P*R**r*K**k*(factorial(r)*factorial(k))**3
 return dict(component=kind,X_order=r,eta_order=k,upper=v,binary=binary,actual_complete_profile_or_scale_certified=False)
def raw_total(recipe,kind,m,binary=False):
 order(m)
 # Validate the stage and component without expanding any large factorial.
 bound(recipe,kind,0,0,True)
 P=recipe['prefactors'][kind];B=max(Q(1),recipe['radial_bases'][kind],recipe['angular_base'])
 if binary:v=lg(max(Q(1),P))+m*lg(B)+(3*m*lg(Q(m)) if m else 0)
 else:
  if m>12:raise ValueError('Use binary raw total-order budget')
  v=P*B**m*factorial(m)**3
 return dict(order=m,upper=v,binary=binary,scope='repair contribution only; raw max over r+k<=m')
def repair_tuple(recipe,h,m,Xhi,Rlo,Rhi):
 order(m)
 if m>11:raise ValueError('Keep high repair input budgets symbolic')
 F=raw_total(recipe,'F',m)['upper'];U=raw_total(recipe,'U',m+1)['upper']
 t=tuple_provider.recipe(h,0,m,Xhi,Rlo,Rhi,F,U)
 return dict(tuple_recipe=t,scope='leading repair contribution only, including its entire axis primitive',actual_complete_leading_tuple_cone_or_scale_certified=False)
def record():
 return dict(status='physical-leading-repair-mixed-Gevrey3-and-radial-obligation',
  Stirling='sum_k |s(r,k)| D^k((k+1)!)^3 <= [8(1+D)]^r(r!)^3, D>=1; r=0 included',
  physical_EU='A*P*[8(1+D)/xmin]^r*K^k*(r!k!)^3; A=162*32600, D=32600',
  physical_F='multiply prefactor by max(1,xmin^-1), radial base becomes [8(1+D)+1]/xmin',
  raw_total='Pkind*max(1,Rkind,K)^m*(m!)^3; U needs m+1 for Q121 primitive',
  pressure='same full five-row target and accepted root remain; a local U plateau does not erase primitive or pressure memory',
  obstruction='angular-only bounds do not control radial jets: a*sin(MX) independent of eta has fixed angular norm but radial derivative a*M',
  next='actual mixed derivatives of compact moving modulation, then all other regions on same hull; bind low-order gates and scales',
  actual_complete_profile_cone_pulse_covariance_residual_or_blowup_certified=False,historical_frozen_records_modified=False)
