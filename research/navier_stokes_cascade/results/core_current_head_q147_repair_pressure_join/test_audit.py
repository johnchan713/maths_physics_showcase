import json,tempfile,unittest
from fractions import Fraction as Q
from pathlib import Path
import q147 as q
import audit

class Tests(unittest.TestCase):
 def post(self):
  p=q.previous.previous.previous.previous
  c=p.previous.selected_caps(2,4,16,Q(1,100),Q(1,2**20),10,100,p.previous.IDENTITY)
  short=p.previous.build(c,Q(1,100),Q(1,2**20),p.previous.IDENTITY)
  args=(Q(1,100),Q(1,2**20),110,4,10,10,p.IDENTITY)
  return p.post_taper(short,p.reference_sources(short,*args),*args)
 def bound(self):
  post=self.post();s=post['state']
  caps=q.previous.cap_tuple(Q(1,8),2,1,Q(1,100),1,1,10,200,16,2010,s['E']['P'],s['invE']['P'],q.previous.IDENTITY)
  return q.previous.bind(post,caps,q.previous.IDENTITY)
 def modulation(self,b):
  r=b['recipe'];p=q.modulation_provider.pullback(2,2,2,1,1,16,q.modulation_provider.IDENTITY)
  return q.modulation_provider.modulation(p,r['field']['F']['P'],r['field']['invF']['P'],r['field']['U']['P'],r['Ay'],r['Be'],q.modulation_provider.IDENTITY)
 def binding(self):
  b=self.bound();return q.bind(b,self.modulation(b),1,2,1,Q(1,4),Q(1,100),q.IDENTITY)
 def test_zero_baseline_is_functional_and_same_gauge(self):
  b=q.baseline(self.bound(),q.IDENTITY);self.assertTrue(b['zero_as_eta_functions'])
  self.assertEqual(b['pairs'],{k:(Q(0),Q(1)) for k in q.ROWS});self.assertIn('Pi_axis',b['pressure'])
 def test_same_q138_memory_has_all_five_rows(self):
  b=self.bound();m=q.memory(b,self.modulation(b),1,2,q.IDENTITY)
  self.assertEqual(tuple(m['rows']),q.ROWS);self.assertTrue(m['same_q138_object']);self.assertTrue(m['full_Pi_and_Pi_eta_memory'])
  self.assertGreater(m['density_prefactors']['Pi'],0);self.assertEqual(m['shifted_derivative_base'],8*m['base'])
 def test_density_envelopes_cover_exact_sample(self):
  b=self.bound();mod=self.modulation(b);m=q.memory(b,mod,1,2,q.IDENTITY);e=mod['prefactors']['E']/9;u=mod['prefactors']['U']-mod['prefactors']['deltaU'];de=mod['prefactors']['deltaE'];du=mod['prefactors']['deltaU']
  exacts=dict(M=du,I=Q(3,2)*de,J=Q(3,2)*((e+de)*(u+du)-e*u),S=(u+du)**2-(e+de)**2/2-(u*u-e*e/2),Pi=((e+de)**2-e*e)/2)
  for row,v in exacts.items():self.assertLessEqual(abs(v),m['density_prefactors'][row])
 def test_signed_q74_and_q132_conventions_agree(self):
  rows=tuple(map(Q,(2,3,5,7,11)));rho=Q(13);e=Q(17);lam=Q(1,100);a=Q(19)
  z=q.terminal_target(rows,rho,e,lam,a);M,I,J,S,P=rows;m=M/(rho*e);i=I/(a*e);j=J/(a*e*e);s=S/(rho*e*e);p=P/(e*e)
  self.assertEqual(z,(-m,(j-m)/lam,-i,s,-p))
 def test_target_and_join_close_without_pressure_reset(self):
  x=self.binding();self.assertTrue(x['same_modulation_memory_certified']);self.assertTrue(x['actual_baseline_defects_certified'])
  c=q.close(x,q.IDENTITY);self.assertTrue(c['same_exact_repair_and_inner_outer_value_cone_certified']);self.assertTrue(c['global_pressure_and_five_moment_identities_certified'])
  self.assertEqual(len(c['global_identities']),5);self.assertIn('integral_X^infty',c['pressure']);self.assertFalse(c['coefficient_reserved_mean_supports_certified'])
 def test_stage_only_or_pressure_dropped_is_rejected(self):
  b=self.bound();mod=self.modulation(b);bad=dict(mod);bad['prefactors']=dict(mod['prefactors']);bad['prefactors']['E']=Q(1,1000)
  with self.assertRaises(ValueError):q.memory(b,bad,1,2,q.IDENTITY)
  x=self.binding();badx=dict(x);badx['memory']=dict(x['memory']);badx['memory']['rows']=q.ROWS[:-1]
  with self.assertRaises(ValueError):q.close(badx,q.IDENTITY)
 def test_00_record(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   out.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+'\n')
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
