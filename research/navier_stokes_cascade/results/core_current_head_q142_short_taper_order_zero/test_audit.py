import json,tempfile,unittest
from fractions import Fraction as Q
from pathlib import Path
import q142 as q
import audit

class Tests(unittest.TestCase):
 def caps(self):return q.selected_caps(2,4,16,Q(1,100),Q(1,2**20),10,100,q.IDENTITY)
 def test_selected_geometry_traces_and_pressure(self):
  c=self.caps();self.assertEqual(c['geometry'],dict(X=2,invX=1,R=4,invR=4))
  self.assertEqual(c['traces']['Pi'],20*16**2+(4/Q(4))*(2*2)**2)
  self.assertEqual(c['trace_cap'],max(c['traces'].values()))
  self.assertTrue(c['full_pressure_trace_included'])
 def test_exponential_caps_and_floor(self):
  c=self.caps();self.assertEqual(c['reference_exp_cap'],729000*8)
  self.assertEqual(c['actual_exp_cap'],3)
  self.assertEqual(c['p1_reference_floor'],Q(1,10))
  self.assertEqual(c['a_floor'],Q(1,10*2**20))
 def test_Q141_map_receives_certified_inputs(self):
  c=self.caps();r=q.build(c,Q(1,100),Q(1,2**20),q.IDENTITY)
  self.assertTrue(r['geometry_trace_exponential_and_a_floor_certified'])
  self.assertTrue(r['recipe']['sources']['all_five_incoming_traces_retained'])
  self.assertTrue(r['recipe']['state']['full_p2_including_Pi_and_Pi_eta'])
  self.assertFalse(r['complete_compact_regional_four_field_maximum_certified'])
 def test_binary_caps_are_not_materialized(self):
  r=q.binary_value_caps(10,100000,100001,q.IDENTITY)
  self.assertEqual(r['base'],2);self.assertGreater(r['exponent'],200000)
  self.assertFalse(r['materialized_fraction'])
 def test_regional_max_min_and_guards(self):
  rs={'short':dict(H=3,Ay=4,Be=5,afloor=Q(1,10),G=6,V0=7),
      'outer':dict(H=8,Ay=3,Be=2,afloor=Q(1,20),G=9,V0=10)}
  z=q.regional_obligation(rs,10)
  self.assertEqual((z['H'],z['afloor'],z['V0']),(8,Q(1,20),10))
  bad={k:dict(v) for k,v in rs.items()};bad['outer']['V0']=11
  with self.assertRaises(ValueError):q.regional_obligation(bad,10)
  with self.assertRaises(ValueError):q.selected_caps(2,4,16,Q(1,100),Q(1,2**19),10,100,q.IDENTITY)
  with self.assertRaises(ValueError):q.selected_caps(2,4,16,Q(1,100),Q(1,2**20),10,100,'wrong')
 def test_00_record(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   out.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+'\n')
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
