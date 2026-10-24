"""Complete initial weighted Haar rows over the symbolic coefficient ring."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
sharp=load('symbolic_haar_sharp','core_sharp_initial_scalar_join/bounds.py')
old=load('symbolic_haar_original','core_initial_haar_row_budget/budgets.py')
EP,order=sharp.EP,sharp.order

def row_budget(n,data):
 order(n)
 if set(data)!=set(old.KEYS):raise ValueError('Every base, fast, slow, covariance and pressure input required')
 b,V,G,F,Fr,u,z,x,w,D,pp,rho,kp,r=(EP(data[k]) for k in old.KEYS)
 d=n+4;vis=2*d*d+d*r+r*r
 S=2*b*x+x*x+Fr*Fr+w;H=b*z+G*x+x*z+Fr*F+w;T=2*V*u+u*u+F*F+w
 gr=d*x+(d+r)*S+d*H+r*T+vis*x
 pm=kp*(gr+rho*pp)
 rt=b*u+V*x+x*u+Fr*F+D;at=G*u+V*z+z*u+F*F+w
 rz=b*z+G*x+x*z+Fr*F+D;az=2*G*z+z*z+F*F+w+pm
 return dict(gr=gr,pm=pm,Etheta=d*u+(d+2*r)*rt+d*at+vis*u,
  Ez=d*z+(d+r)*rz+d*az+vis*z,Srr=S,Srz=H,Stt=T,Rtheta=rt,Atheta=at,Rz=rz,Az=az)

class InitialHaar:
 """Construct row slots from the actual initial norm recipes, not row callbacks.

 higher(N,p) is the remaining weighted Sigma-Sigma0 bound, order3,
 R-chart and pole<=p. c_rho(n) bounds the WHOLE moving c_rho through n.
 """
 def __init__(self,initial,max_output,higher,c_rho):
  if not isinstance(initial,sharp.SharpInitial):raise TypeError('Same-family sharp initial scalar and field join required')
  order(max_output)
  if initial.max_output<max_output:raise ValueError('Sharp scalar/field ledger too short')
  if not all(callable(v) for v in (higher,c_rho)):raise TypeError('Original higher-stress and moving-compensator recipes required')
  self.initial,self.max_output,self.higher,self.c_rho=initial,max_output,higher,c_rho
  self.pole=initial.budget.pole
 def check(self,n):
  order(n)
  if n>self.max_output:raise ValueError('Outside the declared initial Haar output ledger')
 @lru_cache(None)
 def inputs(self,n):
  self.check(n);N=n+3;b=self.initial.budget
  base=b._data(N,2*self.pole);f=b._fast(N);s=self.initial.slow_fields(N,True)
  c=b.primary.covariance_R(N)
  return dict(b=base['b'],V=base['V'],G=base['G'],Ft=f['u'],Fr=f['x'],
   U=s['U'],Zg=s['Z'],X=s['X'],W=EP(c['W']),D=EP(c['Dcurl'])+EP(self.higher(N,self.pole)),
   PiP=self.initial.output(n+1)['P'],Rho=EP(b.fixed(N,self.pole)['rho']),Kp=base['Kp'],r=base['r'])
 @lru_cache(None)
 def rows(self,n):return row_budget(n,self.inputs(n))
 def projected_source(self,n,e,*,pullback,bump):
  """Supply the current row and exact projector moment BEFORE two-edge inversion.

  pullback(n,e): fixed same-hull R -> log-chart jet/multiplier norm.
  bump(n,e): weighted normalized projector bump norm in that chart.
  These fixed inputs do not assert an actual target norm or moment value.
  """
  self.check(n)
  if e not in (1,2):raise ValueError('Original axial or angular row required')
  t=self.initial.output(n+1)
  moment=(n+1)*(t['Jtheta'] if e==2 else t['Jz']+EP(self.c_rho(n+1))*t['P'])
  row=EP(pullback(n,e))*self.rows(n)['Etheta' if e==2 else 'Ez']
  return dict(row=row,moment=moment,H=row+EP(bump(n,e))*moment,
   input_pole=2*self.pole,epsilon_order=2-Q(3,100000),moment_actual_epsilon_order=3,
   signed_primitive_or_root_certified=False)
 def ledger(self):return dict(row_order=2-Q(2,100000),row_pole=2*self.pole,
  field_order_offset=3,pressure_moment_order_offset=1,max_covariance_order=self.max_output+13,
  conditional_initial_Haar_only=True,complete_residual_or_blowup_certified=False)
