"""First REAL signed covariance and complete pressure-change norm substitution."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
signed=load('real_join_signed','core_symbolic_signed_target_join/bounds.py')
generic=load('real_join_pressure','core_generic_signed_covariance_join/budgets.py')
EP,order=signed.EP,signed.order
def rescale(v):
 return EP({k:c.scale(Q(25,16)**k) for k,c in EP(v).c.items()})

class RealCovariance:
 """Only the first phase-locked real signed update; no sourced wave has intervened."""
 def __init__(self,target,max_output):
  if not isinstance(target,signed.SignedTarget):raise TypeError('Same-family original signed-root recipe required')
  order(max_output)
  if target.max_order<max_output+3:raise ValueError('Pressure output n needs signed target through n+3')
  self.target,self.max_output=target,max_output
  self.primary=target.source.haar.initial.budget.primary
  if self.primary.primary.root is not target.root:raise ValueError('The signed and original primary roots must be identical')
 @lru_cache(None)
 def coefficients(self,n):
  order(n)
  if n>self.max_output+2:raise ValueError('Outside the finite covariance/pressure ledger')
  root=self.target.root
  make=lambda f:16*188**n*factorial(n)*rescale(sum((EP(f(j)) for j in range(n+1)),EP()))
  return dict(A=make(root.amplitude),B=make(root.complete_curl),L=make(self.target.signed_amplitude),S=make(self.target.signed_curl))
 @lru_cache(None)
 def covariance(self,n):
  v=self.coefficients(n);A,B,L,S=(v[k] for k in ('A','B','L','S'))
  transfer=EP(self.primary.chart.transfer(n))
  principal=transfer*2**(n+1)*A*L
  error=transfer*2**n*(2*B*S+L*L+S*S)
  q=max(6*n+9,6*n+self.target.pole+6)
  return dict(principal=principal,error=error,complete=principal+error,pole=2*q)
 def pressure(self,n,operators):
  """operators(order,pole) must cover the NEW source and all radial moments.

  Its declared bounds remain lower inequalities, not numerical certificates.
  """
  order(n)
  if n>self.max_output:raise ValueError('Outside the declared pressure output ledger')
  a,b,c=(self.covariance(n+j) for j in range(3));pole=c['pole']
  op=operators(n+1,pole)
  keys={'order','pole','rho','Kp','Kmom','primitive','bump','c_rho'}
  if set(op)!=keys:raise ValueError('Full new-domain pressure, moment, projector and primitive bounds required')
  order(op['order']);order(op['pole'])
  if op['order']<n+1 or op['pole']<pole:raise ValueError('Old pressure proof domain does not cover this covariance update')
  lo=self.primary.chart.Rmin;r=EP(sum(lo**(-j-1) for j in range(n+2)))
  vals={k:EP(op[k]) for k in keys-{'order','pole'}}
  z=generic.pressure_target_budget(n,a['complete'],b['complete'],c['complete'],a['error'],
       r=r,**{k:vals[k] for k in ('rho','Kp','Kmom','primitive','bump')})
  next_change=generic.scalar.wave_change(n+1,b['complete'],c['complete'],vals['Kmom'])
  old=self.target.source.haar.initial
  base=old.output(n);base_next=old.output(n+1)
  changed=dict(P=z['moment_P'],Jtheta=z['Jtheta'],Jz=z['Jz'])
  updated={k:EP(base[k])+changed[k] for k in changed}
  updated_next={k:EP(base_next[k])+next_change[k] for k in changed}
  return dict(changes=z,updated_scalars=updated,
   whole_pressure_change_moment=(n+1)*(next_change['Jz']+vals['c_rho']*next_change['P']),
   whole_updated_pressure_moment=(n+1)*(updated_next['Jz']+vals['c_rho']*updated_next['P']),
   operator_pole=pole,old_Haar_operator_pole=self.target.pole)
 def ledger(self):
  k=Q(1,100000);alpha=2-3*k
  return dict(principal_and_complete_covariance_order=alpha,covariance_error_order=3-6*k,
   Haar_pressure_and_scalar_change_order=alpha,pointwise_pressure_change_order=alpha-k,
   updated_Haar_tangential_order=3-6*k,whole_pressure_moment_order=alpha+1,
   Haar_row_jet_requires_next_target_order=True,
   actual_projector_and_fixed_new_operator_inputs_still_required=True,
   later_complex_or_complete_residual_certified=False)
