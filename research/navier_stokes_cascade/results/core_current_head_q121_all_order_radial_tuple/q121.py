"""All fixed radial tuple jets; actual whole-interval profile bounds remain inputs."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
current=load('q121_current','core_current_head_q120_closed_borel_cutoff/q120.py')
radial=load('q121_radial','core_background_radial_reconstruction/budgets.py')
order=current.order
exact=current.borel.exact

def recipe(h,n,m,Xhi,Rlo,Rhi,Fraw,Uraw):
 order(n);order(m)
 h,Xhi,Rlo,Rhi,Fraw,Uraw=map(exact,(h,Xhi,Rlo,Rhi,Fraw,Uraw))
 if not 0<h<Q(1,100) or not 0<Rlo<=Rhi or not Rlo**2<=2*Xhi<=Rhi**2 or min(Fraw,Uraw)<0:raise ValueError('Same positive hull, h and full raw profile inputs required')
 ell=1-2*h;c=Q(1,2)-h+2*n*h;B=max(Q(1),Xhi)
 tau=min(Q(1,8),Rlo**2/8)
 prefactor=4*Rhi*Fraw+(2+(112+32*c)*B/(ell*Rlo))*Uraw
 return dict(coefficient=n,order=m,tau=tau,prefactor=prefactor,factorial_argument=m,
             geometric_base=1/tau,F_order=m,U_order=m+1,
             whole_axis_interval_required=True,actual_inputs_certified=False,
             infinite_schedule_certified=False,pressure_support_certified=False)

def materialize(r,maximum=40):
 order(maximum)
 if r['order']>maximum:raise ValueError('Keep high factorial as a recipe')
 return r['prefactor']*factorial(r['order'])*r['geometric_base']**r['order']

def tuple_binary_upper(r):
 """No factorial or high power evaluation; zero input gets a padded norm."""
 m=r['order'];lg=current.borel.ceil_log2
 return lg(max(Q(1),r['prefactor']))+(m*lg(Q(m)) if m else 0)+m*lg(r['geometric_base'])

def diagonal_from_profiles(h,n,m,Xhi,Rlo,Rhi,Fraw,Uraw,log2_lower):
 r=recipe(h,n,m,Xhi,Rlo,Rhi,Fraw,Uraw)
 gate=current.deferred_row(h,n,m,tuple_binary_upper(r),log2_lower)
 return dict(tuple_recipe=r,row=gate,actual_profile_and_scale_provenance_certified=False)

def record():
 return dict(status='all-order-axis-memory-radial-tuple-bound',
  tau='min(1,Rlo^2)/8',prefactor='4*Rhi*Fraw+[2+(112+32*(D+2*n*h))*max(1,Xhi)/((1-2*h)*Rlo)]*Uraw',
  raw_tuple_bound='m!*tau^(-m)*prefactor',F_order='m',U_order='m+1 on [0,Xhi]',
  actual_profiles_schedule_pressure_support_or_pulse_certified=False)
