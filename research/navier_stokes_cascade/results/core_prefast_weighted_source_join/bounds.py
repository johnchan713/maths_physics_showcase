"""Explicit primary covariance to complete PREFast zero-Haar source inputs."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
chart=load('prefast_reverse','core_weighted_reverse_chart/bounds.py')
prefix=load('prefast_full_rows','core_signed_prefix_source_budget/budgets.py')
P,order=chart.P,chart.order

class Prefast:
 def __init__(self,primary):
  if not isinstance(primary,chart.ChartedPrimary):raise TypeError('Complete same-hull charted primary required')
  self.primary=primary
 def rows(self,n,data,pole):
  """Remaining data: full Sigma and fixed rho, pressure and whole-moment maps.

  Sigma: epsilon*zeta*delta^-pole. rho: weighted at pole<=pole.
  Kp acts at pole; Kmom bounds the whole radial/Haar moment. No means
  have yet been added. This method is invalid for a postfast prefix.
  """
  order(n);N=n+3
  if set(data)!={'sigma','rho','Kp','Kmom'}:raise ValueError('Full stress and all pressure/moment inputs required')
  if not all(isinstance(v,P) for v in data.values()):raise TypeError('Same-family positive polynomial bounds required')
  if type(pole) is not int or pole<12*N+18:raise ValueError('Common pole must cover the full primary covariance and fixed inputs')
  lo=self.primary.chart.Rmin;B=self.primary.chart.B
  g=P(2+6*B*sum(((j+1)*lo**(-j) for j in range(N+1)),Q(0)))
  r=P(sum((lo**(-j-1) for j in range(N+1)),Q(0)))
  w=self.primary.covariance_R(N)['W'];d=n+4;dr=d*g
  gr=(dr+r)*w+d*w+r*w
  whole=data['Kmom']*gr;pm=data['Kp']*(gr+data['rho']*whole)
  et=(dr+2*r)*(w+data['sigma'])+d*w
  ez=(dr+r)*(w+data['sigma'])+d*(w+pm)
  return dict(gr=gr,P=whole,pm=pm,Etheta=et,Ez=ez,Jtheta=data['Kmom']*w,Jz=data['Kmom']*(w+gr),
              graph=g,r=r,covariance=w,field_order=N,pole=pole,epsilon_order=1-Q(1,100000))
 def fast_inputs(self,m,data,pole):
  """data(N) supplies only the remaining same-hull fixed/background inputs.

  One common pole for BOTH rows; factorials convert normalized row bounds
  to raw sources, and 2 is the zero-Haar projection norm, not a gain.
  """
  order(m)
  a=self.rows(m+4,data(m+7),pole);b=self.rows(m+5,data(m+8),pole)
  return dict(Et_raw=2*factorial(m+4)*a['Etheta'],Ez_raw=2*factorial(m+4)*a['Ez'],
              Ez_radial_raw=2*factorial(m+5)*b['Ez'],covariance_input_order=m+8,pole=pole,
              epsilon_order=a['epsilon_order'],actual_postfast_or_slow_state_certified=False)

def from_original(p,selector_inputs,ell,**kw):
 bundle,b=chart.from_original(p,selector_inputs,ell,**kw)
 return bundle,Prefast(b)
