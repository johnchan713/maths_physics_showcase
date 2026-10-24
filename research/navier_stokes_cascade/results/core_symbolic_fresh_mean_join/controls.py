from bounds import EP,ExpSum,Q,seed,FreshMeans,FIXED,full_rows
def fixture():
 class TracePrimary(seed.chart.ChartedPrimary):
  def __init__(self):
   super().__init__(seed.chart.bridge.WeightedPrimary(seed.chart.bridge.primary.curl.entrance.inverse.roots.fixture()),seed.chart.ReverseChart(1,2));self.seen=[]
  def covariance_R(self,N):self.seen.append(N);return dict(W=seed.P({0:N+1,1:1}))
 primary=TracePrimary();calls=[]
 def fixed(N,p):calls.append(('fixed',N,p));return {k:EP(i+1) for i,k in enumerate(FIXED)}
 def compact(n,p):calls.append(('compact',n,p));return dict(KT=ExpSum([(n+2,0)]),KA=ExpSum([(n+3,0)]),I=6)
 def slow(n):calls.append(('slow',n));return {k:n+i+1 for i,k in enumerate(('H','Kp','Kz','Kg','Kpsi'))}
 b=FreshMeans(primary,fixed,compact,slow,dict(l=2,u=3,Xa=1,Xb=8,t1=primary.primary.root.width),0,162)
 return b,calls
def review():
 checks={};old=seed.prefix
 for n in (0,1,3):
  data={k:old.P({0:i+1,1:n+1}) for i,k in enumerate(old.STATE+old.FIXED+('w',))}
  a=full_rows(n,data);b=old.full_mean(n,data)
  for S in (1,4):
   checks[f'order{n}_S{S}_all_historical_complete_row_outputs']=all(a[k].rational_upper(S)==b[k].evaluate(S) for k in a)
 a=EP({0:ExpSum([(2,3)]),2:ExpSum([(5,7)])});b=EP({1:ExpSum([(3,11)])})
 checks['exact_exponential_polynomial_product']=(a*b).record()=={'1':[['6','14']],'3':[['15','18']]}
 checks['exact_rational_polynomial_copy']=EP(seed.P({0:2,3:4})).rational_upper(2)==34
 budget,calls=fixture();r=budget.rows(0);pre=budget._prefast(0)
 checks['fresh_pressure_and_scalars_differ_from_prefast']=all(r[k].record()!=pre[k].record() for k in ('P','pm','Jtheta','Jz'))
 state=budget._fast(3);d=budget._data(3,324)
 original=old.full_mean(0,{k:old.P({j:v.rational_upper() for j,v in EP(vv).c.items()}) for k,vv in dict(d,**state).items()})
 checks['actual_constructed_fast_state_reproduces_original_rows']=all(r[k].rational_upper(1)==original[k].evaluate(1) for k in r)
 changed=dict(d,**state);changed['Kmom']=0
 checks['whole_moving_pressure_not_deleted']=full_rows(0,changed)['Ez'].record()!=r['Ez'].record()
 original_compact=budget.compact
 def noA(n,p):v=original_compact(n,p);return dict(v,KA=0)
 budget2,_=fixture();budget2.compact=noA
 checks['actual_axial_A1_changes_updated_state']=budget._fast(0)['z'].record()!=budget2._fast(0)['z'].record()
 slow=budget.slow_increment(0);a,b=budget.rows(0),budget.rows(1);k,kn=budget.slow(0),budget.slow(1)
 checks['slow_response_uses_fresh_pressure_and_axial_scalar']=slow['u'].record()==(budget.floor*k['H']*(k['Kp']*a['P']+k['Kz']*a['Jz'])).record()
 checks['radial_slow_response_uses_next_fresh_angular_scalar']=slow['x'].record()==(budget.floor*kn['H']*kn['Kpsi']*b['Jtheta']).record()
 checks['stale_angular_scalar_changes_radial_response']=slow['x'].record()!=(budget.floor*kn['H']*kn['Kpsi']*a['Jtheta']).record()
 final=budget.final(0)
 checks['all_old_fast_means_are_added_not_reset']=all(final[k].record()==(budget._fast(0)[k]+slow[k]).record() for k in final)
 checks['doubled_pressure_operator_pole_requested']=('fixed',3,324) in calls
 checks['safe_covariance_order_twelve_reached']=max(budget.primary.seen)==12
 checks['slow_hull_floor_stays_symbolic']=len(budget.floor.c[0].terms)==1 and next(iter(budget.floor.c[0].terms))>0
 checks['coarse_slow_size_cannot_supply_order_one_uniformly']=all(Q(2)**j==Q(2)**(100000*j*Q(1,100000)) for j in (1,2,4)) and Q(1)-Q(1,100000)<1
 checks['no_later_signed_or_residual_promotion']=budget.ledger()['complete_residual_or_signed_prefix_certified'] is False
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='91b233e043f2bb8f6efd6c81854547e2dc0d5ca0',scope='exact coefficient-ring and manufactured trace diagnostics; no actual continuum supremum certificate')
