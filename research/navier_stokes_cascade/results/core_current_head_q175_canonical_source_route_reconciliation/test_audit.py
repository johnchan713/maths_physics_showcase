import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q175


class Q175Tests(unittest.TestCase):
    def test_statuses(self):
        self.assertTrue(all(q175.statuses().values()))

    def test_polynomial_to_log(self):
        row = q175.polynomial_to_log(7, Q(1, 200), 4096)
        self.assertEqual(row["log_degree"], 14)
        self.assertTrue(row["bound_holds"])
        self.assertLessEqual(row["ratio"], row["uniform_bound"])

    def test_polynomial_to_log_guards(self):
        for args in ((-1, Q(1, 2), 3), (1, 0, 3), (1, Q(1, 2), 0)):
            with self.assertRaises(ValueError):
                q175.polynomial_to_log(*args)

    def test_conditional_route(self):
        row = q175.typed_route(True, True, True, True, False)
        self.assertTrue(row["conditional_typed_source_route_closed"])
        self.assertFalse(row["independently_validated_typed_source_route"])
        self.assertFalse(row["physical_to_canonical_transfer_used"])

    def test_route_requires_each_input(self):
        for missing in range(4):
            flags = [True, True, True, True, False]
            flags[missing] = False
            self.assertFalse(q175.typed_route(*flags)["conditional_typed_source_route_closed"])

    def test_route_guards(self):
        with self.assertRaises(ValueError):
            q175.typed_route(True, True, True, True, 1)

    def test_correction_and_scope(self):
        row = q175.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["correction"]["Q174_same_family_source_frontier_was_stale"])
        self.assertTrue(row["scope"]["conditional_canonical_typed_source_membership_closed_relative_to_named_analytic_providers"])
        self.assertFalse(row["scope"]["provider_theorems_independently_reproved_here"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup_proved"])

    def test_pressure_and_residual_memory(self):
        checks = q175.review()["checks"]
        self.assertTrue(checks["complete_pressure_memory_is_retained"])
        self.assertTrue(checks["complete_finite_residual_is_already_conditionally_joined"])

    def test_invalid_identity(self):
        for fn in (q175.statuses, q175.review):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 420)
        self.assertGreater(record["closure"]["source_bindings"], 5637)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
