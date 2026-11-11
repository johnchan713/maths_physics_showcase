import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q180


class Q180Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q180.statuses().values()))

    def test_exact_cycle_rates(self):
        for j in (0, 1, 19, 1002):
            row = q180.rates(j)
            self.assertEqual(row["B0"], Q(3, 2) - Q(2, 100000))
            self.assertEqual(row["delta"], Q(1, 2) - Q(3, 100000))
            self.assertEqual(row["C_j"], row["B_j"] + Q(1, 2))
            self.assertEqual(row["target_slack"], Q(1, 50000))
        with self.assertRaises(ValueError):
            q180.rates(-1)

    def test_simultaneous_step_closes_source_and_target(self):
        for j in (0, 7, 1002):
            row = q180.coupled_step(j)
            self.assertTrue(row["closes_source_next"])
            self.assertTrue(row["closes_target_next"])
            self.assertTrue(row["closes_coupled_next"])
            self.assertEqual(row["source_minimum_margin"], 0)
            self.assertEqual(row["target_minimum_margin"], Q(1, 50000))

    def test_full_linear_source_rows_are_the_limiter(self):
        row = q180.coupled_step()
        zero_rows = {name for name, value in row["source_margins"].items() if value == 0}
        self.assertEqual(zero_rows, {"sourced_full_linear", "signed_full_linear"})

    def test_complete_memory_and_edge_gate(self):
        row = q180.moment_and_edge_contract()
        for term in ("six covariance", "axial velocity", "Pi1/Pi2", "c_rho P_mean", "zero-Haar", "reserved"):
            self.assertIn(term, row["memory"])
        self.assertIn("zero weighted total", row["primitive_gate"])
        self.assertFalse(row["compactness_alone_suffices"])
        self.assertFalse(row["stage_uniform_constant_claimed"])

    def test_correction_and_scope(self):
        row = q180.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["correction"]["Q160_standalone_sequential_P3_then_P4_reading_is_invalid"])
        self.assertTrue(row["established"]["physical_P3_target_and_both_moments_for_each_finite_cycle"])
        self.assertTrue(row["established"]["physical_P4_generated_source_quotient_and_endpoint_extension"])
        self.assertFalse(row["established"]["constants_uniform_in_cycle_harmonic_jet_or_power"])
        self.assertFalse(row["scope"]["canonical_family_identified"])
        self.assertFalse(row["scope"]["complete_diagonal_endpoint_force_replayed_here"])
        self.assertFalse(row["scope"]["unconditional_complete_residual_force_exclusion_or_blowup"])

    def test_invalid_identity(self):
        for fn in (q180.statuses, q180.moment_and_edge_contract, q180.review):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 425)
        self.assertGreater(record["closure"]["edges"], 1459)
        self.assertGreater(record["closure"]["source_bindings"], 5729)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
