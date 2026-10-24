from types import SimpleNamespace as NS
from bounds import RecomputedSigned,ChosenSource,family_module,moving,real,signed,projector,EP,Q
MovingCompensator=moving.MovingCompensator
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
 c=profile_fixture();family=family_module.RhoSquaredProjectors(c)
 class ManufacturedHaar(projector.haar.InitialHaar):
  # Marked LOWER norm fixture. The new SignedTarget methods are NOT mocked.
  def __init__(self):
   self.initial=c.update.target.source.haar.initial
   self.max_output=3;self.pole=2;self.c_rho=c.c_rho
  def rows(self,n):self.check(n);return dict(Etheta=EP(n+1),Ez=EP(n+2))
 return RecomputedSigned(ManufacturedHaar(),family,0)

def review():
 v=fixture();checks={};c=v.family.compensator
 for n in (0,1,2,3):
  for e in (1,2):
   h=v.source.source(n,e)
   checks[f'order{n}_channel{e}_chosen_log_bump_inserted_once']=v.source.bump(n,e).record()==v.family.logarithmic_bump(n,e).record()
   checks[f'order{n}_channel{e}_source_recomputed_from_current_row_and_moment']=h['H'].record()==(h['row']+v.family.logarithmic_bump(n,e)*h['moment']).record()
   checks[f'order{n}_channel{e}_derived_projector_not_arbitrary_callback']=not h['actual_fixed_bump_seminorms_still_required'] and h['original_fixed_rho_profile_mass_support_and_jets_still_required']
 old=c.update.target
 checks['new_signed_target_not_reused_from_old_projector']=v.target is not old and v.target.source is v.source and v.target.root is old.root
 checks['actual_signed_amplitude_recipe_is_recomputed']=type(v.target) is signed.SignedTarget and v.target.signed_amplitude(0).record()!=old.signed_amplitude(0).record()
 for n in (0,1,2):
  a=v.update.coefficients(n);d=v.covariance(n);tr=EP(v.primary.chart.transfer(n))
  expected=tr*2**n*(2*a['B']*a['S']+a['L']*a['L']+a['S']*a['S'])
  checks[f'order{n}_all_signed_squares_and_curl_cross_terms_retained']=d['error'].record()==expected.record()
 p=v.covariance(2)['pole'];op=v.compensator.operators(v.maps)(1,p)
 checks['new_pressure_uses_native_projectors_not_logarithmic_bumps']=op['bump'].record()==(v.family.weighted(1,1)+v.family.weighted(1,2)).record()
 checks['new_pressure_matches_same_choice_fixed_profile_and_moments']=op['rho'].record()==c.rho(1).record() and op['Kmom']==c.moment(p) and op['c_rho']==c.c_rho(1)
 r=v.pressure(0);a,b,d=(v.covariance(j) for j in range(3))
 direct=real.generic.pressure_target_budget(0,a['complete'],b['complete'],d['complete'],a['error'],r=EP(2),**{k:EP(op[k]) for k in ('rho','Kp','Kmom','primitive','bump')})
 checks['complete_pressure_and_all_scalar_changes_recomputed']=all(r['changes'][k].record()==value.record() for k,value in direct.items())
 previous=v.haar.initial.output(0)
 checks['old_scalar_defects_added_not_reset']=all(r['updated_scalars'][k].record()==(previous[k]+r['changes'][j]).record() for k,j in (('P','moment_P'),('Jtheta','Jtheta'),('Jz','Jz')))
 nxt=real.generic.scalar.wave_change(1,b['complete'],d['complete'],EP(op['Kmom']));previous=v.haar.initial.output(1)
 expected=previous['Jz']+nxt['Jz']+op['c_rho']*(previous['P']+nxt['P'])
 checks['whole_updated_pressure_moment_keeps_original_and_new_terms']=r['whole_updated_pressure_moment'].record()==expected.record()
 checks['old_target_object_still_unchanged']=c.update.target is old and old.source is not v.source
 checks['same_choice_conditional_norms_not_continuum_or_mean_closure']=v.ledger()['new_source_target_curl_covariance_and_pressure_norms_recomputed'] and not v.ledger()['actual_continuum_fields_evaluated'] and not v.ledger()['fresh_fast_slow_or_complete_residual_certified']
 return dict(checks={k:bool(x) for k,x in checks.items()},parent_commit='8076cd797ce1c04b367b1f735764bb728b74d208',scope='same-choice conditional norm recomputation using marked manufactured lower Haar bounds; actual SignedTarget/root/curl and complete pressure recipes are evaluated, not continuum fields')
