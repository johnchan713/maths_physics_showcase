from bounds import P,Q,Prefast,chart,prefix
def fixture():
 b=chart.ChartedPrimary(chart.bridge.WeightedPrimary(chart.bridge.primary.curl.entrance.inverse.roots.fixture()),chart.ReverseChart(1,2))
 return Prefast(b)
def data(N):return {k:P(i+N+1) for i,k in enumerate(('sigma','rho','Kp','Kmom'))}
def review():
 import sympy as s
 checks={};b=fixture();a=b.rows(0,data(3),54)
 inputs={k:prefix.P() for k in prefix.STATE+prefix.FIXED+('w',)}
 vals=dict(data(3),graph=a['graph'],r=a['r'],w=a['covariance'])
 inputs.update({k:prefix.P(dict(v.c)) for k,v in vals.items()})
 original=prefix.full_mean(0,inputs)
 for k in ('gr','P','pm','Etheta','Ez','Jtheta','Jz'):
  checks['historical_zero_mean_'+k]=original[k].c==a[k].c
 changed=data(3);changed['Kmom']=P()
 checks['whole_pressure_moment_cannot_be_deleted']=b.rows(0,changed,54)['Ez'].c!=a['Ez'].c
 t=s.symbols('t',real=True);F=s.cos(t)
 checks['zero_Haar_fluctuation_does_not_vanish']=s.integrate(F,(t,0,2*s.pi))==0 and F.subs(t,0)==1
 eps=s.symbols('eps',positive=True)
 checks['graph_derivative_retains_one_kappa_cost']=1-Q(1,100000)<2-Q(2,100000)
 for N in range(5):
  checks[f'order{N}_graph_and_radius_coefficients']=2+6*4*sum(j+1 for j in range(N+1))==2+2*sum(3*4*s.factorial(j+1)/s.factorial(j) for j in range(N+1))
 # Trace only order routing, not a replacement for the real primary above.
 class TracePrimary(chart.ChartedPrimary):
  def __init__(self):self.chart=chart.ReverseChart(1,2);self.seen=[]
  def covariance_R(self,N):self.seen.append(N);return dict(W=P(N+1))
 tp=TracePrimary();v=Prefast(tp);r=v.fast_inputs(0,data,114)
 aa=v.rows(4,data(7),114);bb=v.rows(5,data(8),114)
 checks['both_covariance_orders_seven_and_eight_requested']=7 in tp.seen and 8 in tp.seen
 checks['raw_tangential_factorial_and_projection_retained']=r['Et_raw'].c==(48*aa['Etheta']).c
 checks['raw_radial_uses_distinct_next_row']=r['Ez_radial_raw'].c==(240*bb['Ez']).c and r['Ez_radial_raw'].c!=(240*aa['Ez']).c
 checks['postfast_claim_rejected']=r['actual_postfast_or_slow_state_certified'] is False
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='0632e406f3244f30273bc04519ccd0b143507bd3',scope='exact manufactured primary/row/routing controls; not actual continuum suprema')
