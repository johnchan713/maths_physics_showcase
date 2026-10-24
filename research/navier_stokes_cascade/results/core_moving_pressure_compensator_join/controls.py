from types import SimpleNamespace as NS
from bounds import MovingCompensator,real,EP,Q,geometry

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
def remaining(n,p):return dict(order=n,pole=p,Kp=3,primitive=5,bump=6)

def review():
 import sympy as s
 from math import comb,factorial
 v=fixture();checks={};t=s.symbols('t')
 for n in range(5):
  ell=sum(v.chart.log_s(j)*t**j for j in range(1,n+1))
  independent=2*s.expand(sum((s.Rational(1,2)*ell)**k/s.factorial(k) for k in range(n+1))).coeff(t,n)
  checks[f'order{n}_dilation_has_no_spurious_log_radial_coordinate']=independent==v.dilation(n)
  inner=sum(v.increment(j)*t**j for j in range(1,n+1))
  outer=sum(v.profile[k]*inner**k for k in range(n+1))
  u=sum(v.dilation(j)*t**j for j in range(n+1))
  independent=sum(s.expand(u*outer).coeff(t,j) for j in range(n+1))
  checks[f'order{n}_full_mixed_compensator_composition']=independent==v.rho_unweighted(n)
 # Exact local jets of the ACTUAL implicit coordinate at s0=1.
 # Polynomial profile is a marked local diagnostic, not a replacement rhohat.
 pc=real.load('moving_pressure_coordinate_checks','core_reserved_geometry_recheck/proof_checks.py')
 Z,T,dR=s.symbols('Z T dR');n=3
 jet=pc.implicit_jet(Q(1,200),Q(1,2),n)
 dj=sum(s.Rational(c.numerator,c.denominator)*Z**i*T**j for (i,j),c in jet.c.items())
 def trunc(f):
  return s.Poly(s.expand(f),Z,T,dR).terms()
 up=s.expand(sum(s.binomial(-s.Rational(1,2),k)*dj**k for k in range(n+1)))
 up=sum(c*Z**a*T**b for (a,b,r),c in trunc(up) if a+b+r<=n)
 x=(2+dR)*up;f=up*(1+x+x*x-x*x*x)
 coeff=trunc(f)
 for m in range(n+1):
  norm=sum(abs(c) for a,c in coeff if sum(a)<=m)
  checks[f'order{m}_actual_implicit_coordinate_mixed_RZT_jets_bounded']=norm<=v.rho_unweighted(m)
 checks['moving_compensator_Z_derivative_is_nonzero']=s.expand(f).coeff(Z,1).subs({T:0,dR:0})!=0
 # Finite endpoint-flat signed diagnostic: mass one does not give L1=1.
 x,R=s.symbols('x R',real=True);phi=(x-2)**4*(3-x)**4
 phi=phi/s.integrate(phi,(x,2,3));hat=phi*(1+10*(x-s.Rational(5,2)))
 mass=s.integrate(hat,(x,2,3));chat=s.integrate(x*x*hat,(x,2,3))/2
 H0=sum(abs(c)*3**a[0] for a,c in s.Poly(hat,x).terms())
 checks['signed_fixed_compensator_has_unit_mass_and_negative_part']=mass==1 and hat.subs(x,s.Rational(9,4))<0
 checks['moment_bound_uses_profile_height_not_unit_mass_as_L1']=abs(chat)<=Q(H0)*9/2 and v.cstar==Q(9,2)*v.profile[0]
 for root in (s.Rational(1,2),s.Rational(1),s.Rational(2)):
  rho=hat.subs(x,R/root)/root
  checks[f'scale{root}_exact_moving_pressure_mass_and_second_moment']=s.integrate(rho,(R,2*root,3*root))==1 and s.integrate(R*R*rho,(R,2*root,3*root))/2==root*root*chat
 for p in (0,1,4,9):
  E=real.signed.maps.edge_absorption(p,v.hull['t1'])
  for j in (-1,0,1,2):
   sup=v.Rmin**-1 if j==-1 else max(Q(1),v.Rmax**j)
   checks[f'pole{p}_radial_weight{j}_covered']=v.moment(p)>=E*(v.Rmax-v.Rmin)*sup
 # Independent integration of a zero-extended C3 polynomial diagnostic.
 # The checks request only jets <=2; all endpoint terms used here vanish.
 F=R*(R-1)**4*(2-R)**4*(1+Z+T*T+R*Z*T)
 for n in range(3):
  source=0
  for i in range(n+1):
   for j in range(n+1-i):
    for k in range(n+1-i-j):
     f=s.diff(F,R,i,Z,j,T,k).subs({Z:0,T:0})/factorial(i)/factorial(j)/factorial(k)
     source+=sum(abs(c)*2**a[0] for a,c in s.Poly(f,R).terms())
  for weight in (-1,0,1,2):
   moment=s.integrate(s.expand(R**weight*F),(R,1,2))
   actual=sum(abs(s.diff(moment,Z,j,T,k).subs({Z:0,T:0}))/factorial(j)/factorial(k) for j in range(n+1) for k in range(n+1-j))
   checks[f'order{n}_weight{weight}_integrated_slow_jets_within_full_source_norm']=actual<=v.moment(0)*source
 pole=v.update.covariance(2)["pole"];op=v.operators(remaining)(1,pole)
 checks['three_independent_operator_inputs_replaced_by_derived_recipes']=op['rho'].record()==v.rho(1).record() and op['Kmom']==v.moment(pole) and op['c_rho']==v.c_rho(1)
 r=v.pressure(0,remaining);direct=v.update.pressure(0,v.operators(remaining))
 checks['complete_pressure_and_old_scalar_additions_preserved']=r['updated_scalars']['P'].record()==direct['updated_scalars']['P'].record() and r['changes']['pm_next'].record()==direct['changes']['pm_next'].record()
 old=v.update.target.source.haar.initial.output(1)
 nextchange=real.generic.scalar.wave_change(1,v.update.covariance(1)['complete'],v.update.covariance(2)['complete'],op['Kmom'])
 expected=EP(old['Jz'])+nextchange['Jz']+v.c_rho(1)*(EP(old['P'])+nextchange['P'])
 checks['whole_current_pressure_product_keeps_old_and_new_defects']=r['whole_updated_pressure_moment'].record()==expected.record()
 checks['interior_floor_uses_moving_support_not_rectangular_union']=v.geom['Rminus']==1 and v.geom['Rplus']==6 and v.floor.record()==EP(real.signed.maps.ExpSum([(1,v.geom['log_weight_comparison'])])).record()
 checks['actual_projector_and_primitive_gaps_remain_explicit']=v.ledger()['pressure_primitive_projector_and_actual_sigma_family_still_required'] and not v.ledger()['full_residual_or_blowup_certified']
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='ab0c64127a9aef62cab16f7ba7d88c97e8b31e33',scope='exact finite implicit/mixed jets and signed compact polynomial diagnostics; original smooth rhohat/support/mass/seminorm premises remain')
