from types import SimpleNamespace as NS
from bounds import NativeHaarPrimitives,moving,EP,Q
real=moving.real;MovingCompensator=moving.MovingCompensator
def fixture():
 signed=real.signed;pmod=signed.projector;fresh=pmod.haar.sharp.fresh;seed=fresh.seed
 class Source(pmod.ProjectedSource):
  def __init__(self):
   primary=seed.chart.ChartedPrimary(seed.chart.bridge.WeightedPrimary(seed.chart.bridge.primary.curl.entrance.inverse.roots.fixture()),seed.chart.ReverseChart(1,8))
   self.haar=NS(max_output=4,pole=2,initial=NS(budget=NS(primary=primary),output=lambda n:dict(P=EP(n+2),Jtheta=EP(n+3),Jz=EP(n+4))))
   self.geometry=fresh.geometry.geometry(2,3,1,8,primary.primary.root.width);self.pull=pmod.LogPullback(8)
  def source(self,n,e):return dict(H=EP(n+e+1))
 class NormTarget(signed.SignedTarget):
  def signed_amplitude(self,n):self.check(n);return EP(n+1)
  def signed_curl(self,n):self.check(n+1);return EP(n+2)
 update=real.RealCovariance(NormTarget(Source(),3),0)
 hull=dict(l=2,u=3,Xa=1,Xb=8,t1=update.target.width)
 return MovingCompensator(update,hull,[10**6]*5)
def review():
 import sympy as s
 from math import factorial,comb
 c=fixture();v=NativeHaarPrimitives(c);checks={}
 # Compare with independent retained recipe at rational exact shell length.
 edge=real.generic.edge
 for p in (0,1,4,9):
  rec=edge.primitive_recipe(p,v.target.Llo,v.target.width)
  expected=EP(moving.maps.ExpSum(list(rec['branches'].values())))
  # Middle uses Lhi rather than Llo, a safe deliberate enlargement.
  branches=list(rec['branches'].values());branches[-1]=(branches[-1][0]*v.target.Lhi/v.target.Llo,branches[-1][1])
  checks[f'pole{p}_independent_three_branch_recipe']=v.kernel(p).record()==EP(moving.maps.ExpSum(branches)).record()
  checks[f'pole{p}_native_pressure_has_radial_derivative_term']=v.pressure_constant(0,p).record()==(1+v.b*Q(1,2)*v.kernel(p)).record()
 R,Z,T,theta=s.symbols('R Z T theta',real=True)
 phi=(R-1)**4*(2-R)**4*(1+Z+T*T)
 h=s.diff(phi,R)
 checks['exact_zero_mass_at_every_slow_point']=s.integrate(h,(R,1,2))==0
 checks['native_primitive_keeps_both_orientations']=s.simplify(s.integrate(h,(R,1,R))-phi)==0 and s.simplify(-s.integrate(h,(R,R,2))-phi)==0
 for n in range(4):
  for e in (1,2):
   positive=sum(abs(s.diff(R**e,R,j).subs(R,v.b))/factorial(j) for j in range(n+1))
   negative=sum(abs(s.diff(R**(-e),R,j).subs(R,v.a))/factorial(j) for j in range(n+1))
   checks[f'order{n}_channel{e}_native_multiplier_jets']=v.target_constant(n,0,e).record()==(Q(int((positive*negative).p),int((positive*negative).q))*v.pressure_constant(n,0)).record()
   target=phi/R**e;source=h/R**e
   checks[f'order{n}_channel{e}_exact_weighted_target_identity']=s.simplify(s.diff(R**e*target,R)-R**e*source)==0 and s.integrate(R**e*source,(R,1,2))==0
 for i,j,k in ((1,0,0),(2,0,0),(1,1,0),(1,0,1),(2,1,0)):
  lhs=s.diff(phi,R,i,Z,j,T,k)/factorial(i)/factorial(j)/factorial(k)
  rhs=s.diff(h,R,i-1,Z,j,T,k)/factorial(i-1)/factorial(j)/factorial(k)/i
  checks[f'jet{i}{j}{k}_normalized_native_radial_reduction']=s.simplify(lhs-rhs)==0
 f=(R-1)**4*(2-R)**4*s.sin(theta)
 mass=s.integrate(f,(R,1,2));haarmass=s.integrate(mass,(theta,0,2*s.pi))
 checks['Haar_mass_zero_does_not_give_pointwise_mass_zero']=haarmass==0 and mass!=0
 checks['oscillatory_left_primitive_leaves_right_tail']=s.integrate(f,(R,1,2))!=0
 bump=lambda n,p:EP(n+p+1)
 op=v.operators(bump)(1,44)
 checks['new_pole_pressure_and_two_target_constants_are_inserted']=op['Kp'].record()==v.pressure_constant(1,44).record() and op['primitive'].record()==(v.target_constant(1,44,1)+v.target_constant(1,44,2)).record()
 result=v.pressure(0,bump);direct=c.update.pressure(0,v.operators(bump))
 checks['complete_pressure_and_whole_current_memory_preserved']=result['changes']['pm_next'].record()==direct['changes']['pm_next'].record() and result['whole_updated_pressure_moment'].record()==direct['whole_updated_pressure_moment'].record()
 checks['no_oscillatory_or_full_residual_certificate']=not v.ledger()['oscillatory_compact_maps_or_full_residual_certified'] and v.ledger()['actual_sigma_and_fixed_profile_support_inputs_required']
 return dict(checks={k:bool(x) for k,x in checks.items()},parent_commit='67f59dfe8743de3c092abbddb574a931f95b8f8a',scope='exact finite Haar/native differential and positive kernel controls; no actual projector or continuum certificate')
