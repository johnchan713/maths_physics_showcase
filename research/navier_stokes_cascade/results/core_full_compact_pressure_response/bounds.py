"""Full auxiliary compact T0/A0 bounds and same-choice covariance response."""
from fractions import Fraction as Q
from math import ceil,comb,factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('compact_pressure_recomputed',HERE.parent/'core_rho_projector_signed_recompute/bounds.py')
recompute=importlib.util.module_from_spec(sp);sp.loader.exec_module(recompute)
EP,order=recompute.EP,recompute.order
exact=recompute.moving.exact

class FullCompactPressure:
 """Lower inputs are ORIGINAL raw profile/chart constants, not fine rows.

 in_chart/out_chart bound raw coefficient jets between R and U charts;
 jacobian bounds dU/dx, U_length the full U integration hull length.
 chi_minus,chi,f0,k0_chi_U are raw through-order fixed-profile bounds.
 k0_chi_U includes the NEXT cutoff derivative. Actual support/gap and
 validity of these lower inequalities remain mathematical premises.
 """
 KEYS={'in_chart','out_chart','jacobian','U_length','chi_minus','chi','f0','k0_chi_U'}
 def __init__(self,state,n,pole,gap,profiles):
  if not isinstance(state,recompute.RecomputedSigned):raise TypeError('Same-choice recomputed signed state required')
  state.check();order(n);order(pole)
  if set(profiles)!=self.KEYS:raise ValueError('Complete original raw T0/A0 profile and chart bounds required')
  self.state,self.order,self.pole=state,n,pole
  self.gap=exact(gap);self.profiles={k:exact(v) for k,v in profiles.items()}
  if min(self.profiles.values())<0 or min(self.profiles[k] for k in ('jacobian','U_length','in_chart','out_chart'))<=0:raise ValueError('Positive hull/chart constants and nonnegative profile bounds required')
  lo,hi=state.target.Llo,state.target.Lhi
  if not 0<self.gap<=lo/2:raise ValueError('Actual plateau gap inside the original shell required')
  r=min(Q(1),lo/4,self.gap);w=state.target.width;a,b=w*w,Q(4);q=max(0,ceil((pole-3)/2))
  def D(c):return sum((Q(comb(q,j)*factorial(j))*(r*r/c)**j for j in range(q+1)),Q(0))/(2*c)
  E=recompute.moving.maps.edge_absorption(pole,w)
  C=EP(recompute.moving.maps.ExpSum([(D(a),b/(lo-r)**2),(D(b),a/(lo-r)**2),(hi*E,(a+b)/(r*r))]))
  Z=EP(recompute.moving.maps.ExpSum([(1,(a+b)/(self.gap*self.gap))]))
  v=self.profiles
  conversion=comb(n+5,5)*factorial(n)*v['in_chart']*v['out_chart']
  self.KT0=conversion*v['jacobian']*2**(3*n)*(v['chi_minus']+v['chi'])*v['f0']*C
  self.KA0=conversion*2**(3*n)*v['U_length']*v['k0_chi_U']*v['f0']*E*Z
 def covariance_response(self,n):
  """Complete CHANGED-wave row response, not the whole current fast source."""
  order(n);s=self.state;s.check();N=n+3
  if self.order<n+1:raise ValueError('Full pressure through n+1 is required')
  c=s.covariance(N)
  if self.pole<c['pole']:raise ValueError('Full compact-pressure pole does not cover the new covariance')
  chart=s.primary.chart;lo=chart.Rmin;B=chart.B
  g=EP(2+6*B*sum((j+1)*lo**(-j) for j in range(N+1)))
  r=EP(sum(lo**(-j-1) for j in range(N+1)));D=c['complete'];d=n+4
  M=EP(s.compensator.moment(self.pole));rho=s.compensator.rho(n+1)
  gr=(d*g+d+2*r)*D;P=M*gr;source=gr+rho*P;pm=self.KT0*source
  return dict(gr=gr,P=P,pm=pm,Etheta=(d*g+d+2*r)*D,
   Ez=(d*g+d+r)*D+d*pm,A0=self.KA0*source,Er=rho*P+self.KA0*source,
   covariance_input_order=N,pressure_order=n+1,
   actual_A0_defect_retained=True,whole_current_fast_source_certified=False)
 def ledger(self):
  return dict(operator_domain='full-auxiliary-shifted-T0-and-A0',raw_profile_order=self.order,
   required_cutoff_profile_order=self.order+1,source_pole=self.pole,
   actual_shift_retained=True,coefficient_derivative_M_loss=0,
   original_raw_profile_chart_and_plateau_inputs_required=True,
   full_current_fast_slow_or_residual_certified=False)
