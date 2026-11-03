"""Real Gamma-tail Gevrey-2 jets; all actual-profile premises stay explicit."""
from fractions import Fraction as Q
from math import factorial
IDENTITY='same-Q76-three-row-relative-heat-root-and-complete-infinite-tail'
def exact(x):
 if isinstance(x,bool) or not isinstance(x,(int,Q)):raise ValueError('Exact rational data required')
 return Q(x)
def order(m):
 if isinstance(m,bool) or not isinstance(m,int) or m<0:raise ValueError('Nonnegative integer order required')
def lg(x):
 x=exact(x)
 if x<=1:return 0
 k=max(0,x.numerator.bit_length()-x.denominator.bit_length())
 return k if x<=2**k else k+1
def parameters(h,Xstar):
 h,Xstar=map(exact,(h,Xstar))
 if not 0<h<=Q(1,100) or Xstar<100:raise ValueError('Original positive heat range and Xstar >=100 required')
 return h,Xstar

def heat(h,X):
 h,X=parameters(h,X)
 return dict(stage='real-Gamma-angular-Gevrey2',prefactor=4*h/X,base=Q(32),
  argument='Z=2*(1-eta^2)/X >=0; eta in [-1,1]',complex_Gamma_extension_claimed=False,
  actual_reference_tail_and_cutoff_certified=False)
def tail(h,Xstar,identity):
 h,Xstar=parameters(h,Xstar)
 if identity!=IDENTITY:raise ValueError('Same reference schedule, complete infinite tail and three-row heat map required')
 g=h**4/Xstar
 return dict(stage='same-heat-full-normalized-target-Gevrey2',g=g,Z=800*g,L=Q(64),
  row_prefactors=dict(Cp=400*h**9/Xstar,S=800*h**9/Xstar,I=40*g),
  unnormalized_rows='16*h*eK^2/XK, 32*h*eK^2, 8*sqrt(2)*eK*sqrt(XK); base 64',
  actual_reference_attenuation_complete_tail_and_target_binding_certified=False)
def root(target,identity):
 if target.get('stage')!='same-heat-full-normalized-target-Gevrey2' or identity!=IDENTITY:raise ValueError('SAME heat target and accepted root required; not the five-row compact repair')
 g,Z,L=map(exact,(target['g'],target['Z'],target['L']))
 if not 0<g<=Q(1,10**10) or Z!=800*g or L!=64:raise ValueError('Original heat target envelope required')
 beta,qnorm=Q(50),Q(100);cstar=6400*g;gamma=2*beta;a=gamma*Z;b=gamma*qnorm
 assert 2*beta*qnorm*cstar<Q(1,2) and 8*a*b<=1
 return dict(stage='same-accepted-three-row-heat-root-Gevrey2',value_cap=cstar,
  positive_prefactor=a,base=8*L,majorant_radius=1/(8*L),majorant_ball=a,
  beta=beta,qnorm=qnorm,gamma=gamma,a=a,b=b,
  actual_original_root_acceptance_and_positive_field_gates_certified=False,
  complex_slow_angle_extension_claimed=False)
def bound(recipe,m,binary=False):
 order(m)
 if recipe.get('stage') not in ('real-Gamma-angular-Gevrey2','same-accepted-three-row-heat-root-Gevrey2','same-relative-heat-patch-swirl-correction-Gevrey2'):raise ValueError('Real heat or same heat root/field recipe required')
 P=recipe['value_cap'] if m==0 and 'value_cap' in recipe else recipe.get('positive_prefactor',recipe.get('prefactor'))
 B=recipe['base']
 if binary:v=lg(max(Q(1),P))+m*lg(max(Q(1),B))+(2*m*lg(Q(m)) if m else 0)
 else:
  if m>12:raise ValueError('Use binary budget at high angular order')
  v=P*B**m*factorial(m)**2
 return dict(order=m,upper=v,factorial_power=2,binary=binary,actual_profile_or_scale_certified=False)
def field(recipe,e_star):
 e_star=exact(e_star)
 if recipe.get('stage')!='same-accepted-three-row-heat-root-Gevrey2' or e_star<=0:raise ValueError('Same relative heat patch and positive scalar entrance required')
 # All three actual unnormalized width-.3 bumps: |beta_j|<=64.
 return dict(stage='same-relative-heat-patch-swirl-correction-Gevrey2',
  prefactor=384*e_star*max(recipe['value_cap'],recipe['positive_prefactor']),base=recipe['base']+2,
  also_bounds_F0_correction_since_sqrt_2X_ge_1=True,
  actual_support_mixed_radial_cone_or_complete_F0_certified=False)
def record():
 return dict(status='real-heat-endpoint-Gevrey2-and-same-compensation-root',
  Gamma='raw d_eta^m(H-1) <= (4h/X)*32^m*(m!)^2, one-sided at poles; no complex Gamma disk',
  multiplier='raw d_eta^m((1+chi*(H-1))^2-1) <= (16h/X)*64^m*(m!)^2',
  infinite_tail='Cp:16h eK^2/XK; S:32h eK^2; I:8sqrt(2)eK sqrt(XK), all times 64^m(m!)^2',
  target='component caps (400h^9,800h^9,40h^4)/Xstar; common Z=800h^4/Xstar, L=64',
  same_root='m>=1 raw <=80000g*512^m*(m!)^2; m=0 retained 6400g; beta50 q100 g=h^4/Xstar',
  relative_patch='raw correction <=30720000 e_star g*514^m*(m!)^2; also F0 since sqrt(2X)>=1',
  endpoint='Gamma Z-Taylor coefficient ratio (m+h)(m+1+h)/(m+1) is unbounded at Z=0',
  pressure='complete tail differences only; common incoming pressure and all other five traces persist locally',
  next='whole exterior F0 angular assembly, then axis/compact cover; Q125 needs F0 through m and U0 through m+2',
  actual_profile_gates_complete_F0_mixed_radial_scale_pulse_residual_or_blowup_certified=False,
  historical_frozen_records_modified=False)
