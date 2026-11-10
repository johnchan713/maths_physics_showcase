"""Selected order-zero caps for Q141's unchanged Q72 short taper."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
previous=load('q142_previous','core_current_head_q141_short_taper_gevrey/q141.py')
exact=previous.exact
IDENTITY='same-selected-Q72-short-taper-and-axis-core-attachment'

def selected_caps(S,Lambda,P,rho,t,Benv,invFcap,identity):
 """Return rational caps proved by Q72 plus the axis/core attachment.

 Benv is the inherited amplitude-independent B.  invFcap is Q72's
 same-field incoming reciprocal cap; it is retained rather than inferred
 from a sampled field.  The selected hierarchy proves t*Benv<1.
 """
 S,Lambda,P,rho,t,Benv,invFcap=map(exact,(S,Lambda,P,rho,t,Benv,invFcap))
 if identity!=IDENTITY:
  raise ValueError('Same selected Q72 taper and axis/core attachment required')
 if S<2 or Lambda<1 or P<16 or not 0<rho<=1 or min(t,Benv,invFcap)<=0:
  raise ValueError('Original positive selected-parameter ranges required')
 if t>Q(1,2**20) or t*Benv>1:
  raise ValueError('Selected C^-20 width and t*B<1 hierarchy required')

 Xa=4/Lambda
 F0=2*S
 U0=2*(9+S/Lambda)
 L0=10*S+2                         # log(10S)+2 <= 10S+2
 h0=1/(1600*S)
 b=rho/(100*S)
 traces=dict(
  M=Xa*U0,
  I=Xa**2*F0,
  J=Xa**2*F0*U0,
  S=Xa*U0**2+Xa**2*F0**2/2,
  Pi=20*P**2+Xa*F0**2)
 T=max(traces.values())

 # C>=2 gives 2t<1/2 and exp(2t)<=1/(1-2t)<2.
 geometry=dict(X=8/Lambda,invX=Lambda/4,R=4,invR=Lambda)
 # Integration by parts costs at most 3 L0_original.  Since e<3,
 # exp(3(log(10S)+2)) < 3^6(10S)^3.
 reference_exp=729000*S**3
 # Ja=-1/2 int kappa p1r, 0<=p1r<B/100 and length 2t.
 # Hence exp(+/- Ja)<=exp(tB/100)<3 (the plus side is actually <=1).
 actual_exp=Q(3)
 p1_floor=2/(5*Lambda)              # p1r>=X/10>=Xa/10
 a_floor=t*p1_floor                 # kappa>=t
 return dict(stage='selected-Q72-short-taper-order-zero-caps',Xa=Xa,
  seed=dict(L0=L0,U0=U0,h0=h0,b=b,F0=F0,invF0=invFcap),
  geometry=geometry,traces=traces,trace_cap=T,angular_base=max(Q(1),1/b),
  reference_exp_cap=reference_exp,actual_exp_cap=actual_exp,
  p1_reference_floor=p1_floor,kappa_floor=t,a_floor=a_floor,
  a_cap=Benv/100,full_pressure_trace_included=True,
  selected_geometry_trace_exponential_and_a_floor_certified=True)

def build(caps,h,t,identity):
 if identity!=IDENTITY or caps.get('stage')!='selected-Q72-short-taper-order-zero-caps':
  raise ValueError('Same selected order-zero caps required')
 h,t=map(exact,(h,t))
 if t!=caps['kappa_floor']:
  raise ValueError('Same selected C^-20 width required')
 s=caps['seed']
 seed=previous.seed(s['L0'],s['U0'],s['h0'],s['b'],s['F0'],s['invF0'],previous.IDENTITY)
 step=previous.switch(t,previous.IDENTITY)
 traces=previous.trace_rows(caps['trace_cap'],caps['angular_base'])
 g=previous.geometry(caps['geometry']['X'],caps['geometry']['invX'],caps['geometry']['R'],caps['geometry']['invR'])
 r=previous.short_taper(seed,step,traces,g,h,t,caps['reference_exp_cap'],caps['actual_exp_cap'],previous.IDENTITY)
 return dict(stage='Q141-short-taper-with-certified-order-zero-inputs',recipe=r,
  a_floor=caps['a_floor'],a_cap=caps['a_cap'],
  geometry_trace_exponential_and_a_floor_certified=True,
  complete_compact_regional_four_field_maximum_certified=False)

def binary_value_caps(Benv,logC_lower,logC_upper,identity):
 """Nonmaterialized powers of two for the inherited reference b/p2 caps."""
 Benv,logC_lower,logC_upper=map(exact,(Benv,logC_lower,logC_upper))
 if identity!=IDENTITY or min(Benv,logC_lower,logC_upper)<=0 or logC_lower>logC_upper or logC_lower<10000*Benv:
  raise ValueError('Same selected C hierarchy and outward logarithmic caps required')
 # |p2r|<C exp(2B); on the short taper Er/Ea<2, so |b|<2C exp(2B).
 # Indeed 2t exp(B)=2 exp(B-20 log C)<1 from log C>=10000B.
 x=1+2*(logC_upper+2*Benv)
 n=(x.numerator+x.denominator-1)//x.denominator
 return dict(stage='binary-order-zero-value-cap',base=2,exponent=n,
  bounds=('abs(b)<2^exponent','abs(reference-p2)<2^exponent'),
  materialized_fraction=False)

def regional_obligation(regions,N):
 N=exact(N)
 if N<1 or N.denominator!=1 or not regions:raise ValueError('Original integer frequency and nonempty compact cover required')
 required={'H','Ay','Be','afloor','G','V0'}
 for name,r in regions.items():
  if set(r)!=required:raise ValueError('Every region needs the same six bound fields')
  H,Ay,Be,afloor,G,V0=map(exact,(r[k] for k in ('H','Ay','Be','afloor','G','V0')))
  if min(H,Ay,Be,afloor)>0 and G>=2 and 0<=V0<=N:continue
  raise ValueError('Regional mixed/floor/root/value/frequency inequality failed')
 return dict(stage='whole-compact-regional-gate',H=max(exact(r['H']) for r in regions.values()),
  Ay=max(exact(r['Ay']) for r in regions.values()),Be=max(exact(r['Be']) for r in regions.values()),
  afloor=min(exact(r['afloor']) for r in regions.values()),G=max(exact(r['G']) for r in regions.values()),
  V0=max(exact(r['V0']) for r in regions.values()),N=N,complete=True)

def record():
 return dict(status='selected-short-taper-order-zero-caps-and-regional-obstruction',
  established='Q72 geometry, five full-memory traces, reference/actual exponential values and a>=2*C^-20/(5*Lambda)',
  correction='the apparent short-taper positive-floor gap is closed by p1r>=X/10 and kappa>=C^-20 from the unchanged axis/core attachment',
  pressure='Pi trace and Pi_eta path remain charged; I and J remain in p1r',
  representation='astronomical value caps are retained as binary logarithmic budgets, with lower/upper log C certificates, rather than materialized integers',
  blocking_inequality='on every Q69-Q71/Q73/Q127 compact region require H,Ay,Be<infinity, afloor>0, G>=max(2,abs(mu),abs(p2),abs(mu*p2)), and V0<=unchanged N; then take max/min on the finite cover',
  next='compress the remaining regional four-field bounds in one convention and apply Q140-Q138 before the full-pressure exact repair/cone join',
  actual_complete_profile_pulse_covariance_residual_or_blowup_certified=False,
  historical_frozen_records_modified=False)
