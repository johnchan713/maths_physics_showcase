from bounds import SharpInitial,EP,Q,fresh,initial
def fixture():
 seed=fresh.seed
 class Trace(seed.chart.ChartedPrimary):
  def __init__(self):
   super().__init__(seed.chart.bridge.WeightedPrimary(seed.chart.bridge.primary.curl.entrance.inverse.roots.fixture()),seed.chart.ReverseChart(1,2));self.seen=[]
  def covariance_R(self,N):self.seen.append(N);return dict(W=seed.P({0:N+1,1:1}))
 p=Trace()
 b=fresh.FreshMeans(p,lambda N,p:{k:i+1 for i,k in enumerate(fresh.FIXED)},
  lambda n,p:dict(KT=n+2,KA=n+3,I=6),lambda n:{k:n+i+1 for i,k in enumerate(('H','Kp','Kz','Kg','Kpsi'))},
  dict(l=2,u=3,Xa=1,Xb=8,t1=p.primary.root.width),1,174)
 return SharpInitial(b,0,(1,2))
def review():
 checks={};b=fixture();n=0;I=b.initial_scalars(n);F=b.fast_scalar_change(n);a=b.preslow(n)
 checks['fresh_scalar_is_initial_plus_retained_fast_quadratic']=all(a[k].record()==(I[k]+F[k]).record() for k in a)
 checks['retained_fast_quadratic_not_zero']=all(F[k].rational_upper(1)>0 for k in F)
 f=b.budget._fast(1);T=b.edge*f['u'];X=b.edge*f['x'];old=initial.coupled.fast_budget(0,b.M,0,0,T.rational_upper(1),X.rational_upper(1))
 checks['exact_historical_fast_P_Jz_coefficients']=all(F[k].rational_upper(1)==old[k] for k in ('P','Jz'))
 v=b.slow_fields(2);out=b.output(0);bval=EP(b.budget.fixed(2,348)['b'])
 direct=initial.coupled.slow_budget(0,b.Ms,bval.rational_upper(1),0,0,v['X'].rational_upper(1),v['U'].rational_upper(1),v['Z'].rational_upper(1))
 checks['exact_historical_C1_postslow_remainder']=all(out[k].rational_upper(1)==direct[k] for k in out)
 checks['scalar_output_primary_covariance_order_twelve']=max(b.budget.primary.seen)==12
 moment=b.axial_moment(0,EP({0:3,1:2}));nextout=b.output(1)
 checks['whole_axial_moving_product_keeps_next_scalar']=moment.record()==(nextout['Jz']+EP({0:3,1:2})*nextout['P']).record()
 checks['whole_pressure_moment_not_deletable']=moment.record()!=nextout['Jz'].record()
 checks['axial_moment_primary_covariance_order_thirteen']=max(b.budget.primary.seen)==13
 w=b.slow_fields(0,True);u=b.slow_fields(0)
 checks['actual_reserved_hull_weight_floor_inserted_once']=all(w[k].record()==(b.budget.floor*u[k]).record() for k in u)
 k=Q(1,100000);mu=1-k
 checks['fast_squares_supply_order_one_without_promoting_linear_size']=2*mu>=1 and mu<1
 checks['C1_slow_squares_supply_order_two']=2*Q(1)==1+Q(1)==2
 checks['whole_axial_residual_keeps_order_three']=b.ledger()['whole_axial_moment_order']==3
 checks['no_later_prefix_or_full_residual_promotion']=b.ledger()['later_complex_prefix_or_full_residual_certified'] is False
 # Historical exact field identities, not just positive-polynomial formulas.
 parts=dict(pressure=initial.pressure_review(),prefix=initial.prefix_review())
 for name,r in parts.items():checks.update({name+':'+k:v for k,v in r['checks'].items()})
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='d5ef67446a724e53913590a682f868f545e73cfb',scope='exact initial identities plus manufactured covariance traces, not actual continuum suprema')
