from types import SimpleNamespace as NS
from bounds import RhoSquaredProjectors,moving,EP,Q
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
 c=fixture();v=RhoSquaredProjectors(c);checks={};x,R,Z,T,S=s.symbols('x R Z T S',real=True)
 phi=(x-2)**4*(3-x)**4;phi=phi/s.integrate(phi,(x,2,3))
 hat=s.expand(phi*(1+10*(x-s.Rational(5,2))));A2=s.integrate(hat*hat,(x,2,3))
 checks['signed_unit_mass_profile_is_not_assumed_positive']=s.integrate(hat,(x,2,3))==1 and hat.subs(x,s.Rational(9,4))<0
 checks['L2_normalization_inverse_has_uniform_support_length_bound']=A2>0 and 1/A2<=1
 c.profile=tuple(Q(sum(abs(co)*comb(a[0],k)*3**(a[0]-k) for a,co in s.Poly(hat,x).terms() if a[0]>=k)) for k in range(5))
 for root in (s.Rational(1,2),s.Rational(1),s.Rational(2)):
  rho=hat.subs(x,R/root)/root
  for e in (1,2):
   sigma=R**(-e)*root*rho*rho/A2
   mass=s.integrate(R**e*sigma,(R,2*root,3*root))
   checks[f'scale{root}_channel{e}_exact_weighted_mass_one']=mass==1
   checks[f'scale{root}_channel{e}_same_moving_support_endpoints']=all(sigma.subs(R,j*root)==0 for j in (2,3))
   F=R+R*R;M=s.integrate(R**e*F,(R,2*root,3*root))
   checks[f'scale{root}_channel{e}_exact_projector_zero_moment_and_idempotence']=s.integrate(R**e*(F-sigma*M),(R,2*root,3*root))==0
 # Exact native mixed jets at an actual implicit-coordinate point, s0=1.
 pc=real.load('rho_squared_coordinate_checks','core_reserved_geometry_recheck/proof_checks.py')
 jet=pc.implicit_jet(Q(1,200),Q(1,2),2)
 for e in (1,2):
  f=R**(-e)*S**(-s.Rational(e+1,2))*hat.subs(x,R/s.sqrt(S))**2/A2
  for n in range(3):
   norm=0
   for i in range(n+1):
    out=pc.Jet(2);power=pc.Jet(2,{(0,0):1})
    for k in range(3):
     val=s.diff(f,R,i,S,k).subs({R:s.Rational(5,2),S:1})/factorial(i)/factorial(k)
     out=out+Q(val)*power;power=power*jet
    norm+=sum(abs(a) for key,a in out.c.items() if sum(key)<=n-i)
   checks[f'order{n}_channel{e}_actual_implicit_native_mixed_jets_bounded']=norm<=v.native(n,e)
  for n in range(4):
   reciprocal=sum(abs(s.diff(R**(-e),R,j).subs(R,c.Rmin))/factorial(j) for j in range(n+1))
   expected=c.hull['u']-c.hull['l'];expected*=Q(reciprocal)*sum(c.dilation(j) for j in range(n+1))*c.rho_unweighted(n)**2
   checks[f'order{n}_channel{e}_full_product_jet_majorant']=v.native(n,e)==expected
   checks[f'order{n}_channel{e}_single_moving_support_floor']=v.weighted(n,e).record()==(c.floor*v.native(n,e)).record()
   pull=real.signed.projector.LogPullback(c.Rmax)
   checks[f'order{n}_channel{e}_single_logarithmic_pullback']=v.logarithmic_bump(n,e).record()==(c.floor*pull.pullback(n,e)*v.native(n,e)).record()
 before=c.update.target
 for e in (1,2):v.logarithmic_bump(2,e)
 checks['inherited_target_is_not_silently_redefined']=c.update.target is before and not v.ledger()['inherited_sigma_identity_proved'] and not v.ledger()['current_signed_target_recomputed']
 checks['nonnegativity_comes_from_square_not_sign_of_rho']=all(q>=0 for q in (hat.subs(x,s.Rational(9,4))**2,hat.subs(x,s.Rational(11,4))**2))
 checks['no_complete_residual_or_blowup_claim']=not v.ledger()['full_residual_or_blowup_certified']
 return dict(checks={k:bool(a) for k,a in checks.items()},parent_commit='b1a956753cf3e96203d151537466b3b153362d42',scope='new admissible projector realization; exact finite signed compact polynomial and implicit-jet diagnostics, not identity with inherited sigma or a continuum field solve')
