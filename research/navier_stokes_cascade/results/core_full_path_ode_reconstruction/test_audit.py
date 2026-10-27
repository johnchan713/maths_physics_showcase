from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
from ode import ambient_lift,picard_factorial_threshold,product_jets,projected_system,residence,review
class FullPathODEReconstructionTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.record=audit.generate()
 def test_projected_and_original_equations(self):
  r=projected_system((2,1,1),(1,-1,0),((0,-2,1),(3,0,-1),(1,2,0)),(3,-2,5),(1,-2,0),7,-2,F(3,5));self.assertEqual(r["constraint_derivative"],0);self.assertTrue(all(x==0 for x in r["residual"]))
 def test_unprojected_normal_source_pressure(self):
  r=projected_system((0,0,1),(0,0,0),((0,0,0),)*3,(0,0,2),(1,0,0),5,1,1);self.assertTrue(all(x==0 for x in r["projected_source"]));self.assertNotEqual(r["pressure_over_i"],0)
 def test_ambient_gaussian_lift(self):
  r=ambient_lift((2,1,1),(1,-1,0),((0,-2,1),(3,0,-1),(1,2,0)),F(3,5),F(4,7));self.assertEqual(r["normal_row"],r["expected"]);self.assertTrue(r["tangent_agreement"])
 def test_picard_does_not_need_smallness(self):
  r=picard_factorial_threshold(700,1100);self.assertLessEqual(r["next_ratio"],F(1,2));self.assertGreater(r["index"],1)
 def test_harmonic_residence(self):
  self.assertEqual(residence(13,1,F(1,7)),13);self.assertEqual(residence(13,3,F(1,7)),F(7,8))
 def test_both_cutoff_products_retain_derivatives(self):
  self.assertNotEqual(product_jets((0,2,-3),(7,11,13))[1],0);self.assertNotEqual(product_jets((0,2,-3),(19,23,29))[1],0)
 def test_scope_corrections(self):
  s=review()["scope"];self.assertTrue(s["canonical_primary_is_homogeneous"]);self.assertTrue(s["sourced_inverse_retains_unprojected_pressure"]);self.assertFalse(s["raw_exit_response_declared_zero"])
 def test_dependencies(self):
  self.assertEqual(len(self.record["dependency_evidence_sha256"]),9);self.assertTrue(self.record["inherited_source_checks"]);self.assertTrue(all(self.record["inherited_source_checks"].values()))
 def test_frozen_record_and_tampering(self):
  frozen=json.loads((audit.HERE/"evidence.json").read_text());audit.validate(frozen);bad=deepcopy(frozen);bad["review"]["scope"]["fixed_projector_contract_closed"]=True
  with self.assertRaises(ValueError):audit.validate(bad)
 def test_no_overwrite(self):
  with tempfile.TemporaryDirectory(dir=audit.HERE) as d:
   p=Path(d)/"kept.json";p.write_text("keep")
   with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
   self.assertEqual(p.read_text(),"keep")
if __name__=="__main__":unittest.main()
