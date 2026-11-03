import unittest,json,math
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import audit,q118
E=q118.EP
text=(audit.HERE.parent/'core_current_head_q117_current_covariant_response/test_audit.py').read_text()
ns=dict(q117=q118.bound,audit=audit,Path=Path,Q=Q,factorial=factorial)
exec(text[text.index('E=q117.EP'):text.index('class Tests(')],ns);old_setup=ns['setup']

def setup(maximum=3):
 old=old_setup();post=old.inputs.source.means.postslow;trace=set()
 b=q118.background.BackgroundBudget(Q(1,200),leading=lambda m:Q(m+2),profiles=lambda n,m:Q(n+m+1),
    cutoff=lambda m:Q(m+1),Rmin=1,Rmax=8,derivative_basis='euler',trace=trace)
 P=q118.background.P
 phase=q118.BackgroundPhase(post,b,maximum=maximum,path=lambda j:P({1:2}) if j==0 else P(3 if j==1 else 0),
    Lmin=1,U=Q(1,2),lambda_max=1,gmin=1)
 w=phase.waves(regular_wave=phase.state.target.root.wave,inverse_normal_square=1)
 response=q118.BackgroundResponse(q118.native.inputs.PathInputs(q118.native.inputs.source.CurrentSource(w)))
 return phase,response,trace

class Tests(unittest.TestCase):
 def test_forward_chart_height_and_source_basis(self):
  phase,_,trace=setup()
  for j in range(4):
   f=q118.charts.forward_recipe(j,8,1)
   actual=factorial(j)*q118.charts.materialize(f)
   self.assertLessEqual(math.log(float(actual)),float(phase.forward_height))
  raw=phase.recipe.normal(0)
  value=phase.normal(0)
  self.assertNotEqual(value.record(),q118.cv(raw).record())
  for name,kind in (('normal','normal'),('matrix','matrix'),('moving_normal','moving_normal'),('profile_rate','profile_rate')):
   norm=sum((Q(5**j,factorial(j))*q118.cv(getattr(phase.recipe,kind)(j)) for j in range(4)),E())
   self.assertEqual(getattr(phase,name)(1).record(),(phase.forward*norm).record())
  self.assertEqual(max(j for name,*rest in trace for j in ([rest[-1]] if name in ('leading','cutoff') else [])),5)

 def test_actual_background_requests_and_nonfrozen_growth(self):
  phase,_,trace=setup();phase.matrix(3)
  self.assertIn(('leading',4),trace);self.assertIn(('cutoff',5),trace)
  self.assertIn(('profile',3,4),trace)
  self.assertTrue(phase.recipe.growth_shape(3).c)
  self.assertTrue(phase.recipe.profile_rate(3).c)
  self.assertEqual(q118.request(0)['background_ceiling'],4)
  self.assertEqual(q118.request(100005)['Borel_Euler_cutoff_ceiling'],100010)
  with self.assertRaises(ValueError):phase.normal(4)

 def test_source_fixed_coefficients_use_current_chart_once_on_return(self):
  phase,_,_=setup();table=phase.linear_table(1,[1,E({1:2})])
  for j,row in enumerate(table):
   N=j+2
   factor=2*q118.cv(phase.state.primary.chart.transfer(N))*16*188**N*factorial(N)
   for key,kind in (('Eik','eikonal'),('bR','bR'),('Bz','Bz')):
    norm=sum((Q(5**k,factorial(k))*phase.jet(kind,k) for k in range(N+1)),E())
    expected=factor*q118.native.inputs.source.waves.rescale(norm)
    self.assertEqual(row[key].record(),expected.record())
  self.assertEqual(table[1]['psi'].record(),E({1:2}).record())

 def test_complete_current_response_and_two_carrier_covariance(self):
  phase,response,trace=setup()
  kw=dict(psi_table=[1,1],tangent=1,length=1,entrance_flat_through=0)
  v=response.build(0,None,[0,1],None,[1,1],**kw)
  self.assertIs(v['source_state'],phase.state)
  self.assertEqual(v['source_ceiling'],1);self.assertEqual(v['background_phase_request']['background_ceiling'],4)
  for j,rate in enumerate(v['derived_rate_jets']):self.assertEqual(rate.record(),response.phase('profile_rate',j).record())
  self.assertEqual(v['complete'][0].record(),(v['principal'][0]+v['curl'][0]).record())
  self.assertFalse(v['phase_provider_inequalities_certified'])
  calls=[];original=response.build
  def call(*a,**kw):calls.append(kw['harmonic']);return original(*a,**kw)
  response.build=call
  factory=q118.bound.covariance.SourcedCovariance(response)
  out=factory.covariance(0,None,[0,1],None,[1,1],**kw)
  self.assertEqual(calls,[1,2]);self.assertEqual(out['covariance_pole'],2*out['wave_pole'])
  self.assertIs(out['source_state'],phase.state);self.assertFalse(out['means_updated'])
  with self.assertRaises(ValueError):factory.moments(0,out)

 def test_bad_domains_or_unbound_inputs_are_rejected(self):
  phase,response,_=setup()
  for linear,rate in (([],None),(None,[])):
   with self.assertRaises(ValueError):response.build(0,linear,[0,1],rate,[1,1],psi_table=[1,1],tangent=1,length=1,entrance_flat_through=0)
  phase.budget.Rmax=7
  with self.assertRaises(ValueError):phase.normal(0)
  phase.budget.Rmax=8;response.state=object()
  with self.assertRaises(ValueError):response.build(0,None,[0,1],None,[1,1],psi_table=[1,1],tangent=1,length=1,entrance_flat_through=0)
  with self.assertRaises(TypeError):q118.BackgroundResponse(None)

 def test_reserve_and_entrance_checks_precede_background_evaluation(self):
  phase,response,trace=setup(maximum=2)
  with self.assertRaises(ValueError):response.build(0,None,[0,1],None,[1,1],psi_table=[1,1],tangent=1,length=1,entrance_flat_through=0)
  self.assertEqual(trace,set())
  phase,response,trace=setup()
  with self.assertRaises(ValueError):response.build(0,None,[0,1],None,[1,1],psi_table=[1,1],tangent=1,length=1,entrance_flat_through=True)
  self.assertEqual(trace,set())

 def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
