from types import SimpleNamespace as NS
from bounds import FullCompactPressure,recompute,EP,Q
RecomputedSigned=recompute.RecomputedSigned;family_module=recompute.family_module;moving=recompute.moving;real=recompute.real;signed=recompute.signed;projector=recompute.projector;MovingCompensator=moving.MovingCompensator
def profile_fixture():
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
def fixture():
 c=profile_fixture();c.profile=c.profile+(c.profile[0],);c.max_order=5;family=family_module.RhoSquaredProjectors(c)
 class ManufacturedHaar(projector.haar.InitialHaar):
  # Marked LOWER norm fixture. The new SignedTarget methods are NOT mocked.
  def __init__(self):
   self.initial=c.update.target.source.haar.initial
   self.max_output=4;self.pole=2;self.c_rho=c.c_rho
  def rows(self,n):self.check(n);return dict(Etheta=EP(n+1),Ez=EP(n+2))
 return RecomputedSigned(ManufacturedHaar(),family,1)

def profiles():return dict(in_chart=3,out_chart=4,jacobian=5,U_length=7,chi_minus=2,chi=3,f0=4,k0_chi_U=5)
def bound_fixture(order=1,pole=None):
 state=fixture();p=state.covariance(3)['pole'] if pole is None else pole
 return FullCompactPressure(state,order,p,Q(1,10),profiles())

def review():
 import sympy as s
 from math import comb,factorial
 checks={};U,V,y,theta,Z,T,a,M=s.symbols('U V y theta Z T a M',real=True)
 wave=(1+s.cos(y))*s.cos(theta)
 cov=s.integrate(wave*wave,(theta,0,2*s.pi))/(2*s.pi)
 haar=s.integrate(cov,(y,0,2*s.pi))/(2*s.pi)
 checks['angular_covariance_is_not_auxiliary_Haar']=s.simplify(cov-(1+s.cos(y))**2/2)==0 and haar==s.Rational(3,4) and s.simplify(cov-haar)!=0
 F=(U+a)**2*(1+Z+T)*s.cos(y+M*a)
 checks['fixed_offset_first_coefficient_derivative_has_no_M_factor']=s.diff(F,U)==2*(U+a)*(1+Z+T)*s.cos(y+M*a)
 checks['fixed_offset_mixed_coefficient_derivatives_have_no_M_factor']=s.diff(F,U,2,Z)==2*s.cos(y+M*a)
 phi=V**4*(1-V)**4;primitive=s.integrate(phi*s.exp(s.I*V),V)
 left=s.exp(s.I*(y-U))*(primitive.subs(V,U)-primitive.subs(V,0))
 right=s.exp(s.I*(y-U))*(primitive.subs(V,1)-primitive.subs(V,U))
 chi=3*U**2-2*U**3;op=(1-chi)*left-chi*right;defect=s.diff(chi,U)*(left+right)
 exact=s.diff(op,U)+s.diff(op,y)-phi.subs(V,U)*s.exp(s.I*y)+defect
 checks['shifted_full_compact_operator_keeps_exact_A0_defect']=s.simplify(s.expand(exact))==0
 lower=s.integrate(phi*(1-V*V/2),(V,0,1)) # cos(V)>=1-V^2/2 on [0,1]
 checks['nonzero_aux_mode_has_zero_Haar_but_nonzero_radial_moment']=s.integrate(s.exp(s.I*y),(y,0,2*s.pi))==0 and lower>0
 checks['A0_cutoff_defect_not_identically_zero']=s.diff(chi,U).subs(U,s.Rational(1,2))!=0 and lower>0
 checks['polynomial_diagnostic_has_zero_extension_jets_at_endpoints']=all(s.diff(phi,V,k).subs(V,end)==0 for end in (0,1) for k in range(4))
 v=bound_fixture();n=v.order;p=v.pole;lo,hi=v.state.target.Llo,v.state.target.Lhi;w=v.state.target.width
 raw=2**(3*n)*(v.profiles['chi_minus']+v.profiles['chi'])*v.profiles['f0']
 karaw=2**(3*n)*v.profiles['U_length']*v.profiles['k0_chi_U']*v.profiles['f0']
 old=moving.maps.weighted_maps(p,lo,w,v.gap,jacobian=v.profiles['jacobian'],KT_profile=raw,KA_profile=karaw)
 r=min(Q(1),lo/4,v.gap);E=moving.maps.edge_absorption(p,w)
 extra=EP(moving.maps.ExpSum([((hi-lo)*E*v.profiles['jacobian']*raw,(w*w+4)/(r*r))]))
 convert=comb(n+5,5)*factorial(n)*v.profiles['in_chart']*v.profiles['out_chart']
 checks['full_T0_matches_independent_original_raw_profile_recipe']=v.KT0.record()==(convert*(EP(old['KT'])+extra)).record()
 checks['full_A0_includes_next_cutoff_profile_and_output_floor']=v.KA0.record()==(convert*EP(old['KA'])).record()
 checks['five_variable_normalized_raw_conversion_retained']=convert==comb(n+5,5)*factorial(n)*3*4
 response=v.covariance_response(0);state=v.state;N=3;chart=state.primary.chart
 fixed=dict(graph=EP(2+6*chart.B*sum((j+1)*chart.Rmin**(-j) for j in range(N+1))),r=EP(sum(chart.Rmin**(-j-1) for j in range(N+1))),rho=state.compensator.rho(1),Kmom=EP(state.compensator.moment(p)),Kp=v.KT0)
 historical=recompute.real.generic.final.fast.covariance_rows(0,state.covariance(N)['complete'],fixed)
 checks['full_changed_covariance_response_matches_original_complete_row_functional']=all(response[k].record()==historical[k].record() for k in ('gr','P','pm','Etheta','Ez'))
 checks['whole_compensated_full_pressure_not_replaced_by_uncompensated_source']=response['pm'].record()!=(v.KT0*response['gr']).record()
 checks['changed_radial_row_keeps_compensator_and_nonzero_A0']=response['Er'].record()==(fixed['rho']*response['P']+v.KA0*(response['gr']+fixed['rho']*response['P'])).record() and bool(response['A0'].record())
 checks['new_covariance_and_pressure_order_are_requested']=response['covariance_input_order']==3 and response['pressure_order']==1
 checks['A0_and_missing_old_full_row_remain_explicit']=response['actual_A0_defect_retained'] and not response['whole_current_fast_source_certified'] and not v.ledger()['full_current_fast_slow_or_residual_certified']
 return dict(checks={k:bool(x) for k,x in checks.items()},parent_commit='c82464277857a178873fa93df08712ad01d727ce',scope='conditional full shifted T0/A0 profile bounds; exact compact polynomial Fourier diagnostics and manufactured lower norm response, not actual continuum raw-profile enclosures')
