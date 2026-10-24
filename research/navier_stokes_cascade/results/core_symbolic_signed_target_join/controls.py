from types import SimpleNamespace as NS
from bounds import SignedTarget,projector,EP,Q,maps
def fixture():
 class Manufactured(projector.ProjectedSource):
  def __init__(self):
   roots=projector.haar.sharp.fresh.seed.chart.bridge.primary.curl.entrance.inverse.roots
   primary=NS(primary=NS(root=roots.fixture()),chart=NS(Rmin=Q(1),Rmax=Q(8)))
   self.haar=NS(max_output=3,pole=2,initial=NS(budget=NS(primary=primary)))
   self.geometry=projector.haar.sharp.fresh.geometry.geometry(2,3,1,8,primary.primary.root.width)
   self.pull=projector.LogPullback(8)
  def source(self,n,e):return dict(H=EP({0:n+e+1,1:1}))
 return SignedTarget(Manufactured(),2)
def review():
 import sympy as s
 checks={};v=fixture()
 checks['same_shell_ratio_recovered_without_new_interval']=v.ratio==8 and v.Llo==Q(14,9) and v.Lhi==7
 old=projector.haar.old.edge
 recipe=old.primitive_recipe(v.pole,v.Llo,v.width)
 checks['symbolic_three_branch_constant_dominates_historical_max']=v.primitive.rational_upper(1)>=old.rational_upper(recipe)
 for n in range(3):
  for e in (1,2):
   expected=v.primitive*v.H(n,e)+(v.H(n-1,e) if n else EP())
   checks[f'order{n}_row{e}_same_two_edge_derivative_recipe']=v.F(n,e).record()==expected.record()
 x=s.symbols('x',positive=True);L=4;a=s.Rational(v.width.numerator,v.width.denominator)**2
 weight_inverse=s.exp(a/x**2+4/(L-x)**2)
 for n in range(5):
  expr=s.cancel(s.diff(weight_inverse,x,n)/weight_inverse)/s.factorial(n)
  for point in (Q(1,2),Q(1),Q(2)):
   val=abs(expr.subs(x,s.Rational(point.numerator,point.denominator)))
   delta=min(Q(1),point,L-point)
   checks[f'order{n}_point{point}_full_reciprocal_weight_jet']=bool(val<=v.reciprocal_weight(n)*delta**(-3*n))
 t=s.symbols('t')
 for n in range(5):
  C=sum(v.increment(j)*t**j for j in range(1,n+1))
  checks[f'order{n}_native_finite_Taylor_coefficients']=all(s.expand(C**k).coeff(t,n)==v.power(n,k) for k in range(n+1))
 checks['nonlinear_native_log_chart_terms_retained']=v.increment(2)>0
 # Independent rational-only signed-root fixture for exact coefficient checks.
 class RationalTarget(SignedTarget):
  def total(self,n):return EP({0:n+1,1:2})
 r=RationalTarget(v.source,2);P=type(r.root.K(0))
 for n in range(3):
  original=r.root.signed(n,lambda j:P({0:j+1,1:2}))
  checks[f'order{n}_original_signed_root_unchanged']=r.signed_amplitude(n).rational_upper(2)==original.evaluate(2)
 for n in (0,1):
  original=sum((r.root.data('curl',j)*(r.root.signed(n-j,lambda k:P({0:k+1,1:2}))+(n-j+1)*r.root.signed(n-j+1,lambda k:P({0:k+1,1:2}))) for j in range(n+1)),P())
  checks[f'order{n}_complete_signed_curl_uses_next_target_jet']=r.signed_curl(n).rational_upper(2)==original.evaluate(2)
 checks['actual_symbolic_target_route_has_finite_coefficients']=bool(v.total(1).record()) and bool(v.signed_amplitude(0).record())
 checks['same_original_root_object_is_retained']=v.root is v.source.haar.initial.budget.primary.primary.root
 checks['full_actual_projector_and_residual_scope_still_open']=v.ledger()['actual_projector_uniform_seminorms_still_required'] and not v.ledger()['actual_signed_wave_full_covariance_or_residual_certified']
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='4fd985d9936c16d06a39437340bf4a8c9cb1b0d3',scope='exact historical recipe, finite reciprocal/Taylor and marked manufactured target controls; no actual continuum projector certificate')
