from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
from gaps import chart_tail,first_exact_dyadic_band,normal_floor,persisted_floors,perturbation_radius,review,selector_margins
class SameBandGapReconstructionTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.record=audit.generate()
 def test_continuum_perturbation_radius(self):
  row=perturbation_radius(F(1,5),3,F(1,7),4,F(1,11),F(1,13),F(1,17),2)
  self.assertGreater(row["radius"],0);self.assertEqual(row["radius"],min(row["rows"].values()))
 def test_persisted_floors(self):
  row=persisted_floors(F(1,5),3,F(1,7),4,F(1,11),F(1,13),F(1,17));self.assertTrue(all(v>0 for v in row.values()))
 def test_exact_dyadic_selector(self):
  ell=first_exact_dyadic_band(F(1,200),(2,17,10**6));self.assertEqual(ell%200,0);self.assertGreater(ell,2)
 def test_one_selector_closes_all_rows(self):
  row=selector_margins(400000,F(1,200),(2,17,10**6),F(1,6057072),1000,F(1,2**1000));self.assertTrue(all(v>0 for v in row.values()))
 def test_chart_and_normal_floor(self):
  chart=chart_tail(3,20);self.assertTrue(chart["inside"] and chart["selected"])
  normal=normal_floor(F(1,32),F(1,101));self.assertEqual(normal["fraction"],F(31,32));self.assertTrue(normal["below_coarse_inverse"])
 def test_flat_edge_correction_and_scope(self):
  scope=review()["scope"];self.assertFalse(scope["stress_magnitude_positive_at_flat_edges"]);self.assertFalse(scope["full_path_ode_closed"])
 def test_dependencies(self):
  self.assertEqual(len(self.record["dependency_evidence_sha256"]),8);self.assertTrue(self.record["inherited_source_checks"]);self.assertTrue(all(self.record["inherited_source_checks"].values()))
 def test_frozen_record_and_tampering(self):
  frozen=json.loads((audit.HERE/"evidence.json").read_text());audit.validate(frozen);bad=deepcopy(frozen);bad["review"]["scope"]["full_path_ode_closed"]=True
  with self.assertRaises(ValueError):audit.validate(bad)
 def test_no_overwrite(self):
  with tempfile.TemporaryDirectory(dir=audit.HERE) as d:
   p=Path(d)/"kept.json";p.write_text("keep")
   with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
   self.assertEqual(p.read_text(),"keep")
if __name__=="__main__":unittest.main()
