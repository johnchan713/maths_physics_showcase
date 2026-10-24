from bounds import EP,Q,LogPullback,ProjectedSource,haar
def fixture():
 sharp=haar.sharp;fresh=sharp.fresh;seed=fresh.seed
 class Trace(seed.chart.ChartedPrimary):
  def __init__(self):super().__init__(seed.chart.bridge.WeightedPrimary(seed.chart.bridge.primary.curl.entrance.inverse.roots.fixture()),seed.chart.ReverseChart(Q(1,4),8))
  def covariance_R(self,N):return dict(W=seed.P(N+1),Dcurl=seed.P(N+2))
 p=Trace();hull=dict(l=2,u=3,Xa=1,Xb=8,t1=p.primary.root.width)
 b=fresh.FreshMeans(p,lambda N,p:{k:i+1 for i,k in enumerate(fresh.FIXED)},lambda n,p:dict(KT=n+2,KA=n+3,I=6),lambda n:{k:n+i+1 for i,k in enumerate(('H','Kp','Kz','Kg','Kpsi'))},hull,1,174)
 h=haar.InitialHaar(sharp.SharpInitial(b,0,(1,2)),0,lambda N,p:EP(N+3),lambda n:EP(n+1))
 return ProjectedSource(h,lambda n,e:EP(n+e+1),hull)
def review():
 import sympy as s
 checks={};t=s.symbols('t');p=LogPullback(2)
 for n in range(5):
  A=sum((int(j==1)+p.log_s(j))*t**j for j in range(1,n+1))
  for b in (Q(1,2),Q(1),Q(3,2)):
   independent=s.expand(sum((b*A)**k/s.factorial(k) for k in range(n+1))).coeff(t,n)
   checks[f'order{n}_parameter{b}_exponential_recurrence']=independent==p.exponential(b,n)
  C=sum(p.increment(j)*t**j for j in range(1,n+1))
  checks[f'order{n}_finite_composition_coefficients']=s.expand(sum(C**k for k in range(n+1))).coeff(t,n)==p.composition(n)
 # Mixed derivatives of a nonlinear actual-form coordinate, with an explicitly
 # tame manufactured scale s=exp(Z+T) at the unit center; both rows.
 x,Z,T,R=s.symbols('x Z T R');F=1+R*Z+R**3*T+Z*T
 Phi=s.exp((x+Z+T)/2)
 def indices(n,dim=3):
  if not dim:yield ();return
  for i in range(n+1):
   for a in indices(n-i,dim-1):yield (i,)+a
 def norm(f,variables,n,point):
  total=0
  for a in indices(n):
   v=f
   for var,k in zip(variables,a):v=s.diff(v,var,k)
   total+=abs(v.subs(point))/s.prod(s.factorial(k) for k in a)
  return total
 for n in range(4):
  K=norm(F,(R,Z,T),n,{R:1,Z:0,T:0})
  for e in (1,2):
   f=Phi**(e+1)*F.subs(R,Phi)/2
   checks[f'order{n}_row{e}_mixed_actual_form_pullback']=norm(f,(x,Z,T),n,{x:0,Z:0,T:0})<=K*p.pullback(n,e).rational_upper(1)
 checks['distinct_axial_and_angular_prefactors']=p.pullback(0,1).rational_upper(1)==2 and p.pullback(0,2).rational_upper(1)==4
 checks['nonlinear_radial_map_second_jet_not_zero']=s.diff(Phi,x,2).subs({x:0,Z:0,T:0})==Q(1,4)
 # Auxiliary normalization-only obstruction, finite-regularity fixture.
 R,Z=s.symbols('R Z',real=True);phi=(R-1)**4*(2-R)**4
 for e in (1,2):
  mass=s.integrate(R**e*phi,(R,1,2));m=s.integrate(R**(e+1)*phi,(R,1,2))/mass
  checks[f'row{e}_positive_family_factor_lower_bound']=1<m<2
  for N in (1,3,8):
   sigma=phi/mass*(1+(R-m)*s.sin(N*Z)/2)
   checks[f'row{e}_frequency{N}_mass_one_with_same_support']=s.integrate(R**e*sigma,(R,1,2))==1
   checks[f'row{e}_frequency{N}_moving_derivative_grows_exactly_N']=s.simplify(s.diff(sigma,Z).subs(Z,0)-N*phi*(R-m)/(2*mass))==0
  checks[f'row{e}_moving_derivative_not_identically_zero']=s.expand(phi*(R-m))!=0
 v=fixture();r=v.source(0,1)
 checks['derived_pullback_substituted_into_actual_Haar_recipe']=r['row'].record()==(v.pull.pullback(0,1)*v.haar.rows(0)['Ez']).record()
 checks['bump_floor_and_native_seminorm_inserted_once']=v.bump(0,1).record()==(v.floor*v.pull.pullback(0,1)*2).record()
 checks['full_pressure_moment_and_scope_retained']=r['moment'].record()==v.haar.projected_source(0,1,pullback=v.pull.pullback,bump=v.bump)['moment'].record() and not r['signed_primitive_or_root_certified'] and r['actual_fixed_bump_seminorms_still_required']
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='e73aac0d9a14ad33eb171db92fbef6c203c576ae',scope='exact formal/mixed-jet and auxiliary-family controls; no actual projector or continuum certification')
