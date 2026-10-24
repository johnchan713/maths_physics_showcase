from types import SimpleNamespace as NS
from bounds import RealCovariance,signed,generic,EP,Q
def fixture():
 pmod=signed.projector;fresh=pmod.haar.sharp.fresh;seed=fresh.seed
 class Source(pmod.ProjectedSource):
  def __init__(self):
   primary=seed.chart.ChartedPrimary(seed.chart.bridge.WeightedPrimary(seed.chart.bridge.primary.curl.entrance.inverse.roots.fixture()),seed.chart.ReverseChart(1,8))
   self.haar=NS(max_output=4,pole=2,initial=NS(budget=NS(primary=primary),output=lambda n:dict(P=EP(n+2),Jtheta=EP(n+3),Jz=EP(n+4))))
   self.geometry=fresh.geometry.geometry(2,3,1,8,primary.primary.root.width);self.pull=pmod.LogPullback(8)
  def source(self,n,e):return dict(H=EP(n+e+1))
 class NormTarget(signed.SignedTarget):
  def signed_amplitude(self,n):self.check(n);return EP({0:n+1,1:1})
  def signed_curl(self,n):self.check(n+1);return EP({0:n+2,1:2})
 return RealCovariance(NormTarget(Source(),3),0)
def operators(n,p):return dict(order=n,pole=p,rho=2,Kp=3,Kmom=4,primitive=5,bump=6,c_rho=7)
def review():
 import sympy as s
 checks={};t,theta=s.symbols('t theta',real=True)
 A=s.Matrix([1+t,2-t,1+t*t]);B=s.Matrix([t,1+t,2-t]);L=s.Matrix([2*t,3-t,1+2*t]);S=s.Matrix([1-t,t*t,2*t])
 U=A*s.cos(theta)+B*s.sin(theta);V=L*s.cos(theta)+S*s.sin(theta)
 delta=s.integrate((U+V)*(U+V).T-U*U.T,(theta,0,2*s.pi))/(2*s.pi)
 principal=(A*L.T+L*A.T)/2;error=(B*S.T+S*B.T+L*L.T+S*S.T)/2
 checks['all_six_real_complete_covariance_entries_exact']=all(s.expand(x)==0 for x in delta-principal-error)
 checks['signed_square_and_curl_terms_not_zero']=all((L*L.T)[i,j]!=0 and (B*S.T+S*B.T)[i,j]!=0 for i in range(3) for j in range(i,3))
 E=s.Matrix([1,2,3]);checks['later_complex_cosine_cross_survives']=any(x!=0 for x in E*L.T+L*E.T)
 k=Q(1,100000);alpha=2-3*k;nu=min(alpha+1-2*k,2*alpha-1,2*alpha-2*k)
 checks['real_signed_square_is_rate_limiter']=nu==3-6*k and nu==2*alpha-1
 checks['scalar_changes_cannot_be_promoted_to_fine_error_order']=alpha<2 and alpha<nu
 for h in (Q(0),Q(1,100)):
  checks[f'common_to_native_target_factor_h{h}']=h*alpha-2*(Q(1,2)+h)==-1-3*h*k
  checks[f'atlas_physical_amplitude_factor_h{h}']=all(0<Q(1,2)+h-h*b<Q(3,4) for b in (Q(1,2),1-k,alpha-Q(1,2),alpha-k))
 v=fixture()
 for n in (0,1,2):
  c=v.coefficients(n);r=v.covariance(n);tr=EP(v.primary.chart.transfer(n))
  expected=tr*2**n*(2*c['B']*c['S']+c['L']*c['L']+c['S']*c['S'])
  checks[f'order{n}_all_retained_error_terms_and_one_chart_transfer']=r['error'].record()==expected.record()
  checks[f'order{n}_complete_change_contains_principal_tensor']=r['complete'].record()==(r['principal']+r['error']).record()
 seen=[]
 def op(n,p):seen.append((n,p));return operators(n,p)
 r=v.pressure(0,op);a,b,c=(v.covariance(j) for j in range(3))
 direct=generic.pressure_target_budget(0,a['complete'].rational_upper(1),b['complete'].rational_upper(1),c['complete'].rational_upper(1),a['error'].rational_upper(1),r=sum(v.primary.chart.Rmin**(-j-1) for j in range(2)),rho=2,Kp=3,Kmom=4,primitive=5,bump=6)
 checks['all_complete_pressure_target_outputs_match_historical_formula']=all(r['changes'][key].rational_upper(1)==val for key,val in direct.items())
 checks['new_pressure_domain_is_requested_explicitly']=seen==[(1,c['pole'])] and c['pole']>v.target.pole
 old=v.target.source.haar.initial.output(0)
 checks['old_scalar_defects_are_added_not_reset']=all(r['updated_scalars'][key].record()==(old[key]+r['changes'][mapkey]).record() for key,mapkey in (('P','moment_P'),('Jtheta','Jtheta'),('Jz','Jz')))
 def noCp(n,p):return dict(operators(n,p),c_rho=0)
 checks['whole_moving_pressure_product_changes_axial_moment']=r['whole_pressure_change_moment'].record()!=v.pressure(0,noCp)['whole_pressure_change_moment'].record()
 checks['updated_whole_pressure_includes_old_defects']=r['whole_updated_pressure_moment'].record()!=r['whole_pressure_change_moment'].record()
 checks['no_later_complex_or_complete_residual_claim']=not v.ledger()['later_complex_or_complete_residual_certified']
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='084a424c17fc348ce6ca98426ac26b71e39e0dc6',scope='exact angular tensors and marked manufactured norm substitutions; not actual continuum projector/operator certificates')
