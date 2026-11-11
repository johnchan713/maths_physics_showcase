import json,tempfile,unittest
from fractions import Fraction as Q
from pathlib import Path
import audit,q173

class Q173Tests(unittest.TestCase):
    def test_statuses(self): self.assertTrue(all(q173.statuses().values()))
    def test_projector_moment_and_support(self):
        row=q173.projector_dilation(Q(1,4),1,2,Q(5,4),Q(7,4),Q(10,7),2)
        self.assertEqual(row["weighted_moment"],1); self.assertTrue(row["strictly_enclosed"])
        self.assertNotEqual(row["actual_support"],row["moving_enclosure"])
    def test_projector_guards(self):
        bad=((1,1,2,Q(5,4),Q(7,4),Q(10,7),3),(2,1,2,Q(5,4),Q(7,4),Q(10,7),1),(1,1,2,1,Q(7,4),Q(10,7),1))
        for args in bad:
            with self.assertRaises(ValueError): q173.projector_dilation(*args)
    def test_pressure_memory(self):
        row=q173.moving_pressure(Q(17,9),Q(1,5),Q(7,3),Q(11,4),Q(13,6))
        self.assertNotEqual(row["crho_Z_P"],0)
        self.assertEqual(row["total_without_epsilon"],row["Jz_Z"]+row["crho_Z_P"]+row["crho_P_Z"])
    def test_bounded_operator_counterexample(self):
        row=q173.bounded_operator_counterexample()
        self.assertEqual(row["uniform_operator_bound"],1)
        self.assertTrue(row["prefix_outputs_strictly_increase"])
        self.assertFalse(row["operator_boundedness_implies_uniform_source_bound"])
    def test_generated_source_gate(self):
        self.assertFalse(q173.generated_source_gate(1,6,Q(1,8))["gate_closed"])
        self.assertTrue(q173.generated_source_gate(1,6,Q(1,8),Q(7,3))["gate_closed"])
    def test_generated_source_guards(self):
        for args in ((-1,0,1),(0,-1,1),(0,0,0)):
            with self.assertRaises(ValueError): q173.generated_source_gate(*args)
    def test_scope(self):
        row=q173.review(); self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q4_fixed_profile_geometry_source_audited_on_Q172_tail"])
        self.assertFalse(row["scope"]["actual_generated_source_weighted_suprema_proved"])
        self.assertIn("sup_{ell>=ell_*}",row["first_open_gate"]["inequality"])
    def test_invalid_identity(self):
        for fn in (q173.statuses,q173.review):
            with self.assertRaises(ValueError): fn(identity="wrong")
    def test_frozen_evidence(self):
        record=audit.generate(); self.assertEqual(record,json.loads((audit.HERE/"evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"],417); self.assertGreater(record["closure"]["source_bindings"],5604)
    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            p=Path(td)/"fresh.json"; p.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+"\n")
            self.assertEqual(json.loads(p.read_text()),audit.generate())
if __name__=="__main__": unittest.main()
