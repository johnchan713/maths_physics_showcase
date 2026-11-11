import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q174


class Q174Tests(unittest.TestCase):
    def test_statuses(self):
        self.assertTrue(all(q174.statuses().values()))

    def test_half_weight_obstruction(self):
        row = q174.half_weight_obstruction(7, 32)
        self.assertGreater(row["lower_log_ratio"], 0)
        self.assertTrue(row["positive_after_four_poles"])
        self.assertFalse(row["blanket_full_weight_follows"])

    def test_half_weight_guards(self):
        for args in ((-1, 3), (1, 0), (1.5, 4)):
            with self.assertRaises(ValueError):
                q174.half_weight_obstruction(*args)

    def test_operator_is_one_way(self):
        pressure = ("Jz_Z", "crho_Z_P", "crho_P_Z", "Pi1", "Pi2", "A1_axial")
        self.assertTrue(q174.operator_gate(True, True, pressure)["operator_output_bound"])
        row = q174.operator_gate(False, True, pressure)
        self.assertFalse(row["operator_output_bound"])
        self.assertFalse(row["operator_bound_implies_source_membership"])

    def test_pressure_memory_required(self):
        row = q174.operator_gate(True, True, ("Jz_Z", "crho_P_Z"))
        self.assertFalse(row["complete_pressure_memory"])
        self.assertFalse(row["operator_output_bound"])

    def test_family_nonidentification(self):
        row = q174.family_join("canonical", "physical")
        self.assertFalse(row["may_transfer_actual_source_bounds"])
        self.assertTrue(q174.family_join("canonical", "canonical")["may_transfer_actual_source_bounds"])
        self.assertTrue(q174.family_join("canonical", "physical", True)["may_transfer_actual_source_bounds"])

    def test_typed_gate(self):
        self.assertFalse(q174.typed_source_gate()["gate_closed"])
        self.assertFalse(q174.typed_source_gate(1, 1, False)["gate_closed"])
        self.assertTrue(q174.typed_source_gate(Q(7, 3), Q(11, 5), True)["gate_closed"])

    def test_typed_gate_guards(self):
        with self.assertRaises(ValueError):
            q174.typed_source_gate(-1, 1, True)
        with self.assertRaises(ValueError):
            q174.typed_source_gate(1, -1, True)

    def test_scope_and_exact_frontier(self):
        row = q174.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q5_operator_mapping_audited_on_Q173_fixed_coefficients"])
        self.assertFalse(row["scope"]["actual_canonical_typed_source_membership_proved"])
        self.assertIn("D^I(f_{j,ell}/P_{j,ell})", row["typed_open_gate"]["harmonic"])
        self.assertIn("zeta^{-1}", row["typed_open_gate"]["mean_target"])

    def test_invalid_identity(self):
        for fn in (q174.statuses, q174.review):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 418)
        self.assertGreater(record["closure"]["source_bindings"], 5623)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
