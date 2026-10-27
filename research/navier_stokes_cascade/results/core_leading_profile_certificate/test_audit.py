from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile, unittest
import audit
from certificate import activation_derivative_cost,bump_supports,collar_radii,fixed_order_strip,reserved_power_law,review,supports_are_ordered_and_reserved

class LeadingProfileCertificateTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.record=audit.generate()
 def test_protected_collar_order(self):
  r=collar_radii(); self.assertLess(r["keep"],r["cut"]); self.assertLess(r["cut"],r["protected_end"])
 def test_fixed_activation_cost(self):
  for order in range(20): self.assertEqual(activation_derivative_cost(order),20*order)
  with self.assertRaises(ValueError): activation_derivative_cost(-1)
 def test_fixed_order_strip_not_uniform_over_order(self):
  rows=[fixed_order_strip(F(1,8),order) for order in range(12)]; self.assertTrue(all(rows[j+1]<rows[j] for j in range(11))); self.assertTrue(all(row>0 for row in rows))
 def test_five_bumps_are_strictly_ordered(self):
  rows=bump_supports(); self.assertEqual(len(rows),5); self.assertTrue(supports_are_ordered_and_reserved()); self.assertEqual(rows[-1][1],F(51,20))
 def test_reserved_pressure_identity(self):
  for lam in (F(1,100),F(1,200),F(1,1000)):
   row=reserved_power_law(lam); self.assertEqual(row["E2_decay"],2*row["E_decay"]); self.assertEqual(row["Pi_y_e2_amplitude"],F(1,2))
 def test_full_pressure_memory_and_heat_scope(self):
  row=review()["analytical_join"]; self.assertTrue(row["full_pressure_memory_retained"]); self.assertFalse(row["global_complex_eta_heat_endpoint_claimed"])
 def test_same_band_contract_remains_next(self):
  scope=review()["scope"]; self.assertTrue(scope["leading_leaf_of_first_background_contract_joined"]); self.assertFalse(scope["same_band_continuum_infima_joined"])
 def test_defining_map_dependencies_are_bound(self):
  self.assertEqual(len(self.record["dependency_evidence_sha256"]),8); self.assertTrue(self.record["inherited_source_checks"]); self.assertTrue(all(self.record["inherited_source_checks"].values()))
 def test_frozen_record_and_tampering(self):
  frozen=json.loads((audit.HERE/"evidence.json").read_text()); audit.validate(frozen); bad=deepcopy(frozen); bad["review"]["analytical_join"]["global_complex_eta_heat_endpoint_claimed"]=True
  with self.assertRaises(ValueError): audit.validate(bad)
 def test_existing_output_is_not_overwritten(self):
  with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
   path=Path(directory)/"kept.json"; path.write_text("keep")
   with self.assertRaises(SystemExit): audit.main(["--output",str(path)])
   self.assertEqual(path.read_text(),"keep")
if __name__=="__main__": unittest.main()
