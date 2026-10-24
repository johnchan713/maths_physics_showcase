import json,unittest
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import audit
from bounds import P,ReverseChart,ChartedPrimary,from_original,bridge
class ReverseChartTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,key):
  v=[v for k,v in self.r['review']['checks'].items() if key in k];self.assertTrue(v,key);self.assertTrue(all(v),key)
 def test_Bell_recurrence_against_independent_expansion(self):self.control('Bell_recurrence')
 def test_fractional_power_jets(self):self.control('normalized_radial_enclosure')
 def test_mixed_flat_weight_raw_norm(self):self.control('mixed_flat_weight_raw_sum')
 def test_factorial_mixed_norm(self):self.control('mixed_factorial_norm')
 def test_all_five_coordinates_retained(self):self.control('five_coordinate')
 def test_unit_transfer_failure(self):self.control('unit_reverse')
 def test_nonlinear_chain_required(self):self.control('nonlinear_reverse')
 def test_raw_factorial_conversion(self):self.control('normalization_factorial')
 def test_complete_row_substitution(self):self.control('complete_initial_row_recipe')
 def test_full_pressure_retained(self):self.control('whole_pressure_moment')
 def test_pole_and_epsilon_unchanged(self):self.control('pole_and_epsilon')
 def test_tensor_transferred_once(self):self.control('tensor_transfer_inserted_once')
 def test_order_zero_is_identity_bound(self):self.assertEqual(ReverseChart(1,2).transfer(0).c,P(1).c)
 def test_hull_validation(self):
  for lo,hi in ((0,1),(-1,1),(2,1)):
   with self.assertRaises(ValueError):ReverseChart(lo,hi)
  with self.assertRaises(TypeError):ReverseChart(0.5,1)
 def test_jet_and_chart_type_validation(self):
  c=ReverseChart(1,2)
  for j in (-1,Q(1,2)):
   with self.assertRaises((ValueError,TypeError)):c.transfer(j)
  with self.assertRaises(ValueError):c.radial(0)
  with self.assertRaises(TypeError):ChartedPrimary(object(),c)
 def test_original_factory_uses_same_operator_hull(self):
  normal=bridge.primary.curl.entrance.reference.normal
  p=normal.frame.inputs(bmin=Q(1,20),bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),cmin=1,cmax=1,Dbox=6,Mz=1)
  sp=normal.selector.fixture();L=normal.normal_length_certificate(p)['length_rule_floor'];sp['length']=Q(-(-L.numerator//L.denominator));ell=normal.selector.select_tail(sp)['ell_tail']
  kw=dict(F=lambda j:P(j+2),G=lambda j:P(j+3),btilde=lambda j:P(j+1),path=lambda j:P({1:1}) if j==0 else P(int(j==1)),tilt=Q(1,100),template=bridge.primary.curl.entrance.inverse.roots.fixture(),profile=lambda j:Q(j+1),dyadic=lambda j:Q(j+1),mesh=lambda j:Q(j+1),transverse=lambda j:Q(j+1),temporal=lambda j:Q(j+1))
  bundle,b=from_original(p,sp,ell,**kw)
  self.assertEqual((b.chart.Rmin,b.chart.Rmax),(bundle[2].Rmin,bundle[2].Rmax));self.assertIs(b.primary.root,bundle[3]);self.assertTrue(all(bundle[4]['gates'].values()))
  self.assertGreater(b.covariance_R(1)['W'].evaluate(1),0)
 def test_all_other_row_inputs_remain_required(self):
  b=ChartedPrimary(bridge.WeightedPrimary(bridge.primary.curl.entrance.inverse.roots.fixture()),ReverseChart(1,2))
  with self.assertRaises(ValueError):b.initial_rows(0,{},P(1))
 def test_frozen_record_and_sources(self):
  audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertEqual(len(self.r['dependency_evidence_sha256']),5)
 def test_scope(self):
  self.assertTrue(all(self.r['gates'].values()));self.assertFalse(self.r['flags']['actual_full_mean_residual_or_blowup_certified']);self.assertTrue(self.r['flags']['other_same_hull_weighted_mean_map_inputs_still_required'])
 def test_tampered_scope_recipe_and_provenance_rejected(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tampered']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
