import json,unittest
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import audit
from bounds import P,primary,rows,WeightedPrimary,from_original
class BridgeTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,key):
  v=[v for k,v in self.r['review']['checks'].items() if key in k];self.assertTrue(v);self.assertTrue(all(v))
 def test_common_physical_weight_coordinates(self):self.control('weight_coordinate_invariant')
 def test_all_six_real_covariances(self):self.control('all_six')
 def test_real_quadrature_and_complex_obstruction(self):self.control('mixed_real');self.control('generic_complex')
 def test_full_curl_fine_order_and_pole(self):self.control('fine_order');self.control('pole')
 def test_Haar_cover_and_Leibniz(self):self.control('normalized_cover');self.control('Leibniz')
 def test_historical_covariance_coefficients(self):self.control('original_weighted_covariance_recipe')
 def test_reverse_chart_gate_not_silently_unit(self):self.control('unit_reverse_chart_transfer');self.control('nontrivial_row_chart_transfer')
 def test_full_mean_row_substitution(self):self.control('complete_mean_row_substitution')
 def test_whole_pressure_moment_retained(self):self.control('whole_pressure_moment_input')
 def test_original_primary_factory_enters_bridge(self):
  normal=primary.curl.entrance.reference.normal
  p=normal.frame.inputs(bmin=Q(1,20),bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),cmin=1,cmax=1,Dbox=6,Mz=1)
  sp=normal.selector.fixture();L=normal.normal_length_certificate(p)['length_rule_floor'];sp['length']=Q(-(-L.numerator//L.denominator));ell=normal.selector.select_tail(sp)['ell_tail']
  kw=dict(F=lambda j:P(j+2),G=lambda j:P(j+3),btilde=lambda j:P(j+1),path=lambda j:P({1:1}) if j==0 else P(int(j==1)),tilt=Q(1,100),template=primary.curl.entrance.inverse.roots.fixture(),profile=lambda j:Q(j+1),dyadic=lambda j:Q(j+1),mesh=lambda j:Q(j+1),transverse=lambda j:Q(j+1),temporal=lambda j:Q(j+1))
  bundle,b=from_original(p,sp,ell,**kw);self.assertIs(b.root,bundle[3]);self.assertTrue(all(bundle[4]['gates'].values()));self.assertGreater(b.covariance(0)['W'].evaluate(1),0)
 def test_missing_mean_and_pressure_inputs_rejected(self):
  b=WeightedPrimary(primary.curl.entrance.inverse.roots.fixture())
  with self.assertRaises(ValueError):b.initial_rows(0,{},P(1),P(1))
  with self.assertRaises(TypeError):WeightedPrimary(object())
  with self.assertRaises(ValueError):b.covariance(-1)
 def test_scope_and_all_gates(self):self.assertTrue(all(self.r['gates'].values()));self.assertFalse(self.r['flags']['actual_full_mean_residual_or_blowup_certified'])
 def test_five_historical_records(self):self.assertEqual(len(self.r['dependency_evidence_sha256']),5)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()))
 def test_tampered_pressure_scope_or_hash_rejected(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tampered']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
