from types import SimpleNamespace as NS
from bounds import CurrentFastSource,prefast_zero_haar_rows,full,fast,EP,Q
recompute=full.recompute;RecomputedSigned=recompute.RecomputedSigned;family_module=recompute.family_module;moving=recompute.moving;real=recompute.real;signed=recompute.signed;projector=recompute.projector;MovingCompensator=moving.MovingCompensator;FullCompactPressure=full.FullCompactPressure
def profile_fixture():
 signed=real.signed;pmod=signed.projector;fresh=pmod.haar.sharp.fresh;seed=fresh.seed
 class Source(pmod.ProjectedSource):
  def __init__(self):
   primary=seed.chart.ChartedPrimary(seed.chart.bridge.WeightedPrimary(seed.chart.bridge.primary.curl.entrance.inverse.roots.fixture()),seed.chart.ReverseChart(1,8))
   primary.primary.root.width=Q(1,2) # moderate manufactured diagnostic width
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
 c=profile_fixture();c.profile=c.profile+(c.profile[0],);c.max_order=5
 family=family_module.RhoSquaredProjectors(c);primary=c.update.primary
 class LowerBudget:
  # Marked LOWER field-norm interface, not an actual initial-field enclosure.
  def __init__(self):self.primary=primary;self.pole=224;self.max_output=5
  def _fast(self,n):return dict(x=EP(n+1),u=EP(n+2),z=EP(3*n+4))
  def fixed(self,n,p):return dict(b=EP(1),V=EP(2),G=EP(3),c=EP(4))
  def slow_increment(self,n):raise AssertionError('Coarse slow path cannot replace SharpInitial.slow_fields')
 class LowerSharp(projector.haar.sharp.SharpInitial):
  def __init__(self):self.budget=LowerBudget();self.max_output=4
  def output(self,n):return dict(P=EP(n+2),Jtheta=EP(n+3),Jz=EP(n+4))
  def slow_fields(self,n,weighted=False):
   if weighted is not True:raise AssertionError('The actual moving-hull weighted slow fields are required')
   return dict(X=EP(n+5),U=EP(n+6),Z=EP(n+7))
 class LowerHaar(projector.haar.InitialHaar):
  def __init__(self):self.initial=LowerSharp();self.max_output=4;self.pole=224;self.c_rho=c.c_rho
  def rows(self,n):self.check(n);return dict(Etheta=EP(n+1),Ez=EP(n+2))
 state=RecomputedSigned(LowerHaar(),family,1)
 p=state.covariance(3)['pole']
 profiles=dict(in_chart=3,out_chart=4,jacobian=5,U_length=7,chi_minus=2,chi=3,f0=4,k0_chi_U=5)
 return CurrentFastSource(FullCompactPressure(state,1,p,Q(1,10),profiles),Q(1,12))
def equal(a,b):
 # Compare exact coefficients without serializing potentially huge integers.
 return {k:v.terms for k,v in EP(a).c.items()}=={k:v.terms for k,v in EP(b).c.items()}
def review():
 import sympy as s
 checks={};y=s.symbols('y',real=True)
 def haar(f):return s.integrate(f,(y,0,2*s.pi))/(2*s.pi)
 def Qcirc(f):return s.expand(f-haar(f))
 pre=[2+s.sin(y),3+s.cos(y)];tail=s.sin(2*y)
 rf=[1+s.cos(y),2+s.sin(y)];rs=[3+s.sin(2*y),4+s.cos(2*y)];wc=[5+s.cos(3*y),6+s.sin(3*y)]
 time=[-Qcirc(pre[0]),-Qcirc(pre[1])+tail]
 final=[pre[i]+time[i]+rf[i]+rs[i]+wc[i] for i in (0,1)]
 expected=[Qcirc(rf[0]+rs[0]+wc[0]),tail+Qcirc(rf[1]+rs[1]+wc[1])]
 checks['exact_current_Qcirc_cancellation_identity']=all(s.simplify(Qcirc(final[i])-expected[i])==0 for i in (0,1))
 checks['original_A1_tail_is_axial_and_zero_Haar']=haar(tail)==0 and Qcirc(final[0])==expected[0]
 checks['deleting_original_axial_tail_breaks_identity']=s.simplify(Qcirc(final[1])-Qcirc(rf[1]+rs[1]+wc[1]))==tail
 checks['prefast_Haar_terms_are_annihilated_not_bounded_as_fine_inputs']=haar(pre[0])==2 and haar(pre[1])==3 and all(haar(Qcirc(f))==0 for f in final)
 p0=prefast_zero_haar_rows(0,3,4,5,6)
 checks['original_prefast_zero_Haar_source_program_is_explicit']=equal(p0['Etheta'],90) and equal(p0['Ez'],2235) and p0['Qcirc_norm_factor']==2 and p0['covariance_input_order']==3
 v=fixture();out=v.rows(0);s0=v.state;f=v.budget._fast(3);sl=v.initial.slow_fields(3,True)
 checks['initial_fast_axial_field_is_not_silently_identified_with_angular']=not equal(f['z'],f['u'])
 checks['sharp_weighted_slow_recipe_is_used']=all(equal(out['original_slow_difference'][k],fast.difference(0,f,dict(x=sl['X'],u=sl['U'],z=sl['Z']),dict(b=EP(1),V=EP(2),G=EP(3),c=EP(4),graph=EP(2+6*s0.primary.chart.B*10),r=EP(4),rho=s0.compensator.rho(1),Kp=v.pressure.KT0,Kmom=EP(s0.compensator.moment(v.pressure.pole))))[k]) for k in ('Etheta','Ez','pm'))
 for k in ('Etheta','Ez'):
  checks[k+'_regular_source_retains_both_differences_and_complete_covariance']=equal(out['regular'][k],2*(out['original_fast_difference'][k]+out['original_slow_difference'][k]+out['complete_signed_covariance_difference'][k]))
 checks['full_pressure_and_viscosity_differences_remain_nonzero']=all(out[g][k].c for g in ('original_fast_difference','original_slow_difference') for k in ('pressure_gradient','viscosity','pm'))
 t=out['axial_A1_tail']
 checks['original_tail_request_and_strict_gain_are_exact']=t['integration_order']==99999 and t['raw_source_order']==200001 and t['epsilon_order']=='99999/50000'
 checks['original_primary_source_and_full_pressure_requests_not_truncated']=t['source_program']['covariance_order']==200004 and t['source_program']['full_pressure_order']==200002 and t['minimum_primary_source_pole']==2400066
 checks['small_existing_pole_does_not_cover_original_tail_source']=not t['existing_pole_covers_this_source'] and not t['actual_tail_polynomial_materialized']
 checks['axial_cutoff_gap_not_copied_from_pressure_gap']=t['axial_cutoff_gap']==Q(1,12) and v.pressure.gap==Q(1,10)
 checks['original_tail_program_not_generic_later_stage_callback']='ORIGINAL primary prefast Ez bound at J' in t['original_factors'] and not t['source_program']['independent_fine_source_callback']
 r=v.rates();checks['coarse_A1_tail_prevents_unconditional_fine_rate']=r['whole_source_without_raw_A1_gain']<r['regular_current_source'] and r['whole_source_if_original_A1_obligations_hold']==r['regular_current_source'] and not r['actual_fine_source_or_mean_closure_certified']
 req=v.fast_recipe(0);checks['next_fast_new_layer_requests_and_actual_axial_repair_retained']=req['maximum_new_signed_covariance_order']==8 and req['maximum_new_signed_target_order']==9 and req['minimum_recomputed_pressure_max_output']==6 and req['sufficient_profile_interface_order']==10 and 'gd-A1 gd' in req['actual_new_axial_velocity']
 checks['old_nested_primary_requests_not_mislabeled_as_new_covariance']=req['original_primary_from_initial_fields_order']==18 and req['original_primary_from_A1_tail_order']==200009
 checks['complete_fine_axial_source_not_falsely_certified']=not out['actual_whole_fine_source_certified'] and out['axial_source_regular_plus_unmaterialized_A1']
 return dict(checks={k:bool(x) for k,x in checks.items()},parent_commit='035c285108ddf1b42c764ab8aa5e7f053c36a916',scope='exact Qcirc identity, explicit original-tail obligation and manufactured lower field norms at moderate diagnostic width; no actual continuum fields or fine-source certificate')
