from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from controls import fixture
from bounds import NativeHaarPrimitives
class PrimitiveTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_two_edge_kernel(self):self.control('three_branch');self.control('radial_derivative_term')
 def test_native_target_multipliers(self):self.control('multiplier_jets');self.control('weighted_target_identity')
 def test_native_primitive_orientation_and_slow_mass(self):self.control('every_slow');self.control('both_orientations')
 def test_normalized_radial_jets(self):self.control('radial_reduction')
 def test_oscillatory_counterexample(self):self.control('pointwise_mass');self.control('right_tail')
 def test_pressure_and_memory_join(self):self.control('are_inserted');self.control('current_memory')
 def test_scope(self):self.control('no_oscillatory')
 def test_domains(self):
  v=NativeHaarPrimitives(fixture())
  with self.assertRaises(ValueError):v.kernel(-1)
  with self.assertRaises(ValueError):v.pressure_constant(5,4)
  with self.assertRaises(ValueError):v.target_constant(0,4,3)
  with self.assertRaises(TypeError):NativeHaarPrimitives(object())
  with self.assertRaises(TypeError):v.operators(None)
  c=fixture();c.hull['Xb']*=2
  with self.assertRaises(ValueError):NativeHaarPrimitives(c)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tamper(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
