"""Mixed fixed-zeta primitive bounds from six mixed scalar jets, under original contracts."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
previous=load('q139_previous','core_current_head_q138_mixed_modulation_pullback/q138.py')
scalar=load('q139_scalar','core_current_head_q130_mixed_loop_binding/q130.py')
exact=previous.exact;order=previous.order;lg=previous.lg
IDENTITY='same-six-mixed-scalars-full-pressure-augmented-v-and-fixed-origin-weighted-means'
def phase_disk(disk,identity):
 if identity!=IDENTITY or disk.get('stage')!='analytic-scalar-parameter-disk':raise ValueError('Same Q130 scalar/angle disk and original mean/phase contracts required')
 r,W,F,G=map(exact,(disk['r'],disk['W'],disk['F'],disk['G']))
 if not 0<r<=Q(1,4) or min(W,F)<1 or G<2:raise ValueError('Original disk value bounds required')
 d=1/(8*G**2);t=min(r/4,d*r/(64*W));s=min(r/4,d*t*r/(768*W));delta=d*t/8
 lip=(12*W*s+2*W*t)/(r*d);forcing=(96*W*s/r+delta)/d
 assert lip<Q(1,16) and forcing<=t/4 and forcing+lip*t<t/2
 return dict(stage='same-joint-scalar-phase-inverse-disk',scalar_radius=s,phase_radius=delta,theta_ball=t,density_floor=d,
  primitive_complex_cap=3*F,F=F,r=r,W=W,lipschitz=lip,forcing=forcing,
  off_constraint_complex_periodicity_assumed=False,actual_real_caps_floors_mean_and_disk_binding_certified=False)
def bind(H,Ay,Be,disk,identity):
 H,Ay,Be=map(exact,(H,Ay,Be))
 if identity!=IDENTITY or disk.get('stage')!='same-joint-scalar-phase-inverse-disk' or min(H,Ay,Be)<1:raise ValueError('Actual SIX-scalar mixed jets and SAME joint inverse/mean contracts required')
 s,delta=map(exact,(disk['scalar_radius'],disk['phase_radius']));L=8*(1+H/s)
 return dict(stage='same-six-scalar-bound-mixed-zero-mean-primitives',M=2*disk['primitive_complex_cap'],A=L*Ay,B=L*Be,C=1/delta,
  scalar_jet_prefactor=H,scalar_y_base=Ay,scalar_eta_base=Be,
  actual_mixed_scalar_state_root_augmented_variance_and_pressure_binding_certified=False,
  complex_slow_profile_extension_claimed=False)
def pullback(recipe,V0,N,identity):
 if identity!=IDENTITY or recipe.get('stage')!='same-six-scalar-bound-mixed-zero-mean-primitives':raise ValueError('Same SIX-scalar bound primitives required')
 return previous.pullback(recipe['M'],recipe['A'],recipe['B'],recipe['C'],V0,N,previous.IDENTITY)
def bound(recipe,i,k,j,binary=False):
 for n in (i,k,j):order(n)
 if recipe.get('stage')!='same-six-scalar-bound-mixed-zero-mean-primitives':raise ValueError('Same fixed-zeta mixed primitive recipe required')
 M,A,B,C=(recipe[n] for n in ('M','A','B','C'))
 if binary:v=lg(max(Q(1),M))+i*lg(A)+k*lg(B)+j*lg(C)+(3*i*lg(Q(i)) if i else 0)+(3*k*lg(Q(k)) if k else 0)+(j*lg(Q(j)) if j else 0)
 else:
  if max(i,k,j)>12:raise ValueError('Use binary mixed primitive budget')
  v=M*A**i*B**k*C**j*(factorial(i)*factorial(k))**3*factorial(j)
 return dict(y_order=i,eta_order=k,phase_order=j,upper=v,binary=binary,actual_profile_cone_or_scale_certified=False)
def record():
 return dict(status='six-mixed-scalar-jets-to-fixed-origin-mixed-primitives',
  joint='d=1/(8G^2); t=min(r/4,dr/(64W)); s=min(r/4,dtr/(768W)); delta=d*t/8',
  inverse='fixed-origin phase, complex scalar/phase perturbations; Lipschitz<1/16, image<t/2',
  primitive='fixed integrals over [0,1] define analytic extension of original primitive with cap3F; no off-constraint periodicity premise',
  mixed='M=6F, A=8Ay(1+H/s), B=8Be(1+H/s), C=delta^-1; raw <=M A^i B^k(i!k!)^3 j!C^j',
  scalar_input='sum of raw six-scalar jets/(i!k!) <= H Ay^i Be^k(i!k!)^2; all genuine mixed orders required',
  pressure='p is full p2 with incoming Pi/Pi_eta; v is augmented variance; original weighted means retained on actual real profile',
  next='bind actual mixed six-scalar jets from PREMOD state/cutoff/signed-root, then apply Q138 and join whole hull',
  actual_complete_profile_pulse_covariance_residual_or_blowup_certified=False,historical_frozen_records_modified=False)
