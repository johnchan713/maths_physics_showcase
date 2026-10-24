from bounds import InitialHaar,row_budget,old,sharp,EP,Q
def fixture():
 fresh=sharp.fresh;seed=fresh.seed
 class Trace(seed.chart.ChartedPrimary):
  def __init__(self):
   super().__init__(seed.chart.bridge.WeightedPrimary(seed.chart.bridge.primary.curl.entrance.inverse.roots.fixture()),seed.chart.ReverseChart(1,2));self.seen=[]
  def covariance_R(self,N):self.seen.append(N);return dict(W=seed.P({0:N+1,1:1}),Dcurl=seed.P({0:N+2,1:2}))
 p=Trace();b=fresh.FreshMeans(p,lambda N,p:{k:i+1 for i,k in enumerate(fresh.FIXED)},
  lambda n,p:dict(KT=n+2,KA=n+3,I=6),lambda n:{k:n+i+1 for i,k in enumerate(('H','Kp','Kz','Kg','Kpsi'))},
  dict(l=2,u=3,Xa=1,Xb=8,t1=p.primary.root.width),1,174)
 initial=sharp.SharpInitial(b,0,(1,2))
 return InitialHaar(initial,0,lambda N,p:EP(N+3),lambda n:EP({0:n+1,1:1}))
def review():
 checks={}
 for n in (0,1,3):
  data={k:old.P({0:i+1,1:n+1}) for i,k in enumerate(old.KEYS)}
  a=row_budget(n,data);b=old.row_budget(n,data)
  for S in (1,4):checks[f'order{n}_S{S}_all_complete_original_row_polynomials']=all(a[k].rational_upper(S)==b[k].evaluate(S) for k in a)
 v=fixture();data=v.inputs(0);N=3;b=v.initial.budget
 checks['all_original_input_slots_filled']=set(data)==set(old.KEYS)
 checks['sharp_slow_fields_substituted_on_actual_weighted_hull']=all(data[k].record()==v.initial.slow_fields(N,True)[j].record() for k,j in (('U','U'),('Zg','Z'),('X','X')))
 checks['complete_fast_state_and_actual_A1_kept']=data['Ft'].record()==b._fast(N)['u'].record() and data['Fr'].record()==b._fast(N)['x'].record()
 c=b.primary.covariance_R(N)
 checks['selected_curl_square_inserted_once_plus_higher_stress']=data['D'].record()==(EP(c['Dcurl'])+v.higher(N,v.pole)).record()
 checks['pressure_defect_uses_order_n_plus_one']=data['PiP'].record()==v.initial.output(1)['P'].record()
 checks['max_primary_covariance_order_is_n_plus_thirteen']=max(b.primary.seen)==13
 r=v.rows(0);changed=dict(data,PiP=0);without=row_budget(0,changed)
 checks['whole_pressure_defect_changes_axial_row']=r['Ez'].record()!=without['Ez'].record()
 checks['pressure_change_does_not_spuriously_change_angular_row']=r['Etheta'].record()==without['Etheta'].record()
 kw=dict(pullback=lambda n,e:EP(1),bump=lambda n,e:EP(n+e+1))
 axial=v.projected_source(0,1,**kw);angular=v.projected_source(0,2,**kw)
 t=v.initial.output(1)
 checks['axial_projector_keeps_full_moving_pressure_product']=axial['moment'].record()==(t['Jz']+v.c_rho(1)*t['P']).record()
 checks['angular_projector_has_its_own_distinct_scalar']=angular['moment'].record()==t['Jtheta'].record() and angular['moment'].record()!=axial['moment'].record()
 checks['weighted_projected_source_is_row_plus_bump_moment']=axial['H'].record()==(axial['row']+2*axial['moment']).record()
 checks['fixed_pullback_is_not_silently_unit']=v.projected_source(0,1,pullback=lambda n,e:EP(2),bump=kw['bump'])['row'].record()!=(axial['row']).record()
 checks['no_signed_root_or_complete_residual_promotion']=not axial['signed_primitive_or_root_certified'] and not v.ledger()['complete_residual_or_blowup_certified']
 checks['signed_original_order_and_moment_margin']=axial['epsilon_order']==2-Q(3,100000)<v.ledger()['row_order']<3
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='2701c3c23d605e065a0b6b308d475b8e9420f8c4',scope='exact historical row and manufactured norm-trace controls; not actual continuum suprema')
