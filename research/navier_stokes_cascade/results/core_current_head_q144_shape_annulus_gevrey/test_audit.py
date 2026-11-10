import json,tempfile,unittest
from fractions import Fraction as Q
from pathlib import Path
import q144 as q
import audit

class Tests(unittest.TestCase):
 def post(self,P=16):
  p=q.previous;c=p.previous.selected_caps(2,4,P,Q(1,100),Q(1,2**20),10,100,p.previous.IDENTITY)
  short=p.previous.build(c,Q(1,100),Q(1,2**20),p.previous.IDENTITY)
  args=(Q(1,100),Q(1,2**20),110,4,10,10,p.IDENTITY)
  src=p.reference_sources(short,*args)
  return p.post_taper(short,src,*args)
 def recipe(self,P=16):
  post=self.post(P);s=post['state']
  Ecap=max(Q(100),2*s['E']['P']);inv=max(Q(100),2*s['invE']['P'])
  return q.shape_annulus(post,Q(1,8),2,1,Q(1,100),10,Ecap,inv,200,max(16,P),2010,q.IDENTITY)
 def test_exact_real_floor(self):
  r=q.floor_certificate(1000,1)
  self.assertEqual(r['shape'],Q(12493,15625))
  worst=q.floor_certificate(1,Q(1,100))
  self.assertEqual(worst['shape'],Q(1298,3125))
  self.assertEqual(worst['floor'],Q(1298,3125))
  self.assertGreater(r['earlier_annulus'],Q(7,10))
  self.assertFalse(r['high_order_complex_root_cap_used_for_floor'])
 def test_mixed_same_root_and_complete_memory(self):
  r=self.recipe();self.assertEqual(r['reference_root']['real_value_cap'],q.annular.D)
  self.assertTrue(r['all_five_actual_incoming_rows_reintegrated'])
  self.assertTrue(r['full_p2_Pi_and_Pi_eta_retained'])
  self.assertFalse(r['field']['order_dependent_norm_exponentiated'])
  self.assertGreaterEqual(r['Ay'],1);self.assertGreaterEqual(r['Be'],1)
 def test_pressure_and_endpoint_memory_response(self):
  a=self.recipe(16);b=self.recipe(32)
  self.assertGreater(b['state']['p']['P'],a['state']['p']['P'])
  self.assertGreater(b['incoming_actual_moments']['Pi']['P'],a['incoming_actual_moments']['Pi']['P'])
 def test_precise_gate_request(self):
  g=q.gate_request(self.recipe());self.assertFalse(g['complete'])
  self.assertGreater(g['afloor'],0);self.assertEqual(len(g['missing']),4)
  self.assertIn('Eexp>=exp(G)',g['missing'])
 def test_guards(self):
  with self.assertRaises(ValueError):q.floor_certificate(1,Q(1,50))
  post=self.post();s=post['state'];E=2*s['E']['P'];I=2*s['invE']['P']
  with self.assertRaises(ValueError):q.shape_annulus(post,Q(1,8),2,1,Q(1,100),10,E,I,200,16,2009,q.IDENTITY)
  with self.assertRaises(ValueError):q.shape_annulus(post,Q(1,8),2,1,Q(1,100),10,s['E']['P']/2,I,200,16,2010,q.IDENTITY)
  with self.assertRaises(ValueError):q.shape_annulus(post,Q(1,8),2,1,Q(1,100),10,E,I,200,16,2010,'wrong')
 def test_00_record(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   out.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+'\n')
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
