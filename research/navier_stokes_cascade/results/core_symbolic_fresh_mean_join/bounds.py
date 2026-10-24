"""Positive exponential-coefficient polynomials and fresh primary mean join."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
seed=load('fresh_join_seed','core_prefast_weighted_source_join/bounds.py')
maps=load('fresh_join_maps','core_symbolic_weighted_maps/bounds.py')
geometry=load('fresh_join_geometry','core_reserved_geometry_recheck/proof_checks.py')
order,exact,ExpSum=seed.order,maps.exact,maps.ExpSum

class EP:
 """Finite positive polynomials in S with exact c*exp(e) coefficients."""
 def __init__(self,v=0):
  if isinstance(v,EP):self.c=dict(v.c);return
  if isinstance(v,ExpSum):self.c={0:v} if v.terms else {};return
  if hasattr(v,'c'):v=v.c # exact copies from independent historical P classes
  if not isinstance(v,dict):v={0:v}
  self.c={}
  for k,c in v.items():
   order(k);c=c if isinstance(c,ExpSum) else ExpSum([(c,0)])
   if c.terms:self.c[k]=c
 def __add__(self,v):
  v=EP(v);out=dict(self.c)
  for k,c in v.c.items():out[k]=out.get(k,ExpSum([]))+c
  return EP(out)
 __radd__=__add__
 def __mul__(self,v):
  v=EP(v);out={}
  for i,a in self.c.items():
   for j,b in v.c.items():out[i+j]=out.get(i+j,ExpSum([]))+a*b
  return EP(out)
 __rmul__=__mul__
 def record(self):return {str(k):v.record() for k,v in sorted(self.c.items())}
 def rational_upper(self,S,max_power=10000):
  S=exact(S)
  if S<0:raise ValueError('Nonnegative scale required')
  return sum((c.rational_upper(max_power)*S**k for k,c in self.c.items()),Q(0))

def full_rows(n,data):
 """Same positive historical full_mean formula in the extended coefficient ring."""
 order(n)
 if set(data)!=set(seed.prefix.STATE+seed.prefix.FIXED+('w',)):raise ValueError('All complete state/background/pressure entries required')
 x,u,z,b,V,G,sigma,g,c,r,rho,kp,km,w=(EP(data[k]) for k in seed.prefix.STATE+seed.prefix.FIXED+('w',))
 d=n+4;dr=g*d;time=d*(1+c);vis=d*d*(g*g+1)+r*dr+r*r
 srr=2*b*x+x*x+w;srz=b*z+G*x+x*z+w;stt=2*V*u+u*u+w
 gr=time*x+(dr+r)*srr+d*srz+r*stt+vis*x
 defect=km*gr;pm=kp*(gr+rho*defect)
 rt=b*u+V*x+x*u+w+sigma;rz=b*z+G*x+x*z+w+sigma
 at=G*u+V*z+z*u+w;az=2*G*z+z*z+w+pm
 et=time*u+(dr+2*r)*rt+d*at+vis*u
 ez=time*z+(dr+r)*rz+d*az+vis*z
 return dict(gr=gr,P=defect,pm=pm,Etheta=et,Ez=ez,
             Jtheta=km*(G*u+V*z+z*u+w),Jz=km*(2*G*z+z*z+w+gr))

FIXED=('b','V','G','sigma','c','rho','Kp','Kmom')
class FreshMeans:
 """One primary -> fast -> FRESH scalars -> slow, with no signed wave yet.

 fixed(N,p) supplies original background jets and fixed pressure/moment
 maps at pole p. compact(m,p) supplies actual weighted KT,KA and I,
 I>=c_i^-1/S. slow(n) supplies H,Kp,Kz,Kg,Kpsi for total field jets.
 There are no independent current-state, source-row or scalar callbacks.
 """
 def __init__(self,primary,fixed,compact,slow,hull,max_output,pole):
  if not isinstance(primary,seed.chart.ChartedPrimary):raise TypeError('Same-family charted complete primary required')
  order(max_output)
  if type(pole) is not int or pole<12*(max_output+12)+18:raise ValueError('Common pole must cover all requested primary orders')
  if not all(callable(f) for f in (fixed,compact,slow)):raise TypeError('Original fixed input recipes required')
  if set(hull)!={'l','u','Xa','Xb','t1'}:raise ValueError('Whole original moving potential hull and width required')
  if exact(hull['t1'])!=primary.primary.root.width:raise ValueError('Retain the SAME original primary flat width')
  self.primary,self.fixed,self.compact,self.slow=primary,fixed,compact,slow
  self.max_output,self.pole=max_output,pole
  self.hull=geometry.geometry(**hull)
  self.floor=EP(ExpSum([(1,self.hull['log_weight_comparison'])]))
 def _data(self,N,p):
  d=self.fixed(N,p)
  if set(d)!=set(FIXED):raise ValueError('Full original base and pressure/moment norms required')
  lo=self.primary.chart.Rmin;B=self.primary.chart.B
  d={k:EP(v) for k,v in d.items()}
  d.update(graph=EP(2+6*B*sum((j+1)*lo**(-j) for j in range(N+1))),
           r=EP(sum(lo**(-j-1) for j in range(N+1))),w=EP(self.primary.covariance_R(N)['W']))
  return d
 @lru_cache(None)
 def _prefast(self,n):
  return full_rows(n,dict(self._data(n+3,self.pole),x=0,u=0,z=0))
 @lru_cache(None)
 def _fast(self,m):
  a=self._prefast(m+4);b=self._prefast(m+5)
  ka=self.compact(m,self.pole);kt=self.compact(m+1,self.pole)
  if set(ka)!=set(('KT','KA','I')) or set(kt)!=set(('KT','KA','I')):raise ValueError('Both weighted compact operators and time floor required')
  c=comb(m+5,5)*24*2**(m+4)*EP(ka['I'])*EP({1:1})
  tang=2*factorial(m+4)*c*(a['Etheta']+(1+EP(ka['KA']))*a['Ez'])
  rad=2*factorial(m+5)*comb(m+5,5)*24*2**(m+5)*EP(ka['I'])*EP({1:1})*EP(kt['KT'])*b['Ez']
  return dict(x=rad,u=tang,z=tang)
 @lru_cache(None)
 def rows(self,n):
  order(n)
  if n>self.max_output+1:raise ValueError('Outside the declared finite output ledger')
  return full_rows(n,dict(self._data(n+3,2*self.pole),**self._fast(n+3)))
 def slow_increment(self,n):
  order(n)
  if n>self.max_output:raise ValueError('Outside the declared finite output ledger')
  a,b=self.rows(n),self.rows(n+1);k,kn=self.slow(n),self.slow(n+1)
  if any(set(v)!=set(('H','Kp','Kz','Kg','Kpsi')) for v in (k,kn)):raise ValueError('Every original slow response/profile norm required')
  k={key:EP(v) for key,v in k.items()};kn={key:EP(v) for key,v in kn.items()}
  return dict(u=self.floor*k['H']*(k['Kp']*a['P']+k['Kz']*a['Jz']),
              z=self.floor*k['H']*k['Kg']*a['Jtheta'],
              x=self.floor*(n+1)*kn['H']*kn['Kpsi']*b['Jtheta'])
 def ledger(self):return dict(prefast_pole=self.pole,postfast_pole=2*self.pole,slow_increment_pole=0,
                             max_covariance_order=self.max_output+12,epsilon_mean_order=Q(99999,100000),
                             complete_residual_or_signed_prefix_certified=False)
 def final(self,n):
  inc=self.slow_increment(n);old=self._fast(n)
  return {k:old[k]+inc[k] for k in ('x','u','z')}
