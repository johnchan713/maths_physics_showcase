import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q196


class Q196Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q196.statuses().values()))

    def test_stage_rates(self):
        row = q196.rates()
        self.assertEqual(row["B_initial_source_log"], Q(99999, 100000))
        self.assertEqual(row["B_initial_source_pure_after_kappa"], Q(49999, 50000))
        self.assertEqual(row["D_after_first_sourced_correction"], Q(74999, 50000))
        self.assertLess(row["B_initial_source_log"], row["D_after_first_sourced_correction"])

    def test_linear_ledger(self):
        row = q196.linear_ledger()
        self.assertEqual(row["minimum_gain"], Q(49999, 100000))
        self.assertEqual(row["source_order"], Q(99999, 100000))
        self.assertTrue(row["pressure_gradient_retained"])
        self.assertTrue(row["mixed_viscosity_retained"])
        self.assertTrue(row["source_normal_pressure_retained"])

    def test_all_initial_source_groups(self):
        row = q196.source_groups()
        self.assertTrue(row["all_equal_B"])
        self.assertEqual(set(row["orders"].values()), {Q(99999, 100000)})
        self.assertEqual(row["nonzero_harmonics"], (1, -1, 2, -2))

    def test_polynomial_absorption(self):
        row = q196.absorption(1)
        self.assertEqual(row["degree_ceiling"], 1)
        self.assertEqual(row["uniform_upper"], 3200000000000000)
        self.assertFalse(row["exact_exponent_supremum_follows_without_reserve"])
        self.assertTrue(q196.absorption(0)["exact_exponent_supremum_follows_without_reserve"])

    def test_quotient_and_endpoint_contracts(self):
        row = q196.quotient_contract(5, Q(3, 2))
        self.assertEqual(row["pure_exponent"], Q(49999, 50000))
        self.assertTrue(row["same_domain_raw_source_jets"])
        self.assertTrue(row["smooth_mixed_endpoint_zero_extension"])
        self.assertFalse(row["source_divided_by_cutoff_or_amplitude"])
        endpoint = q196.endpoint_contract()
        self.assertTrue(endpoint["source_over_P_mixed_endpoint_jets_vanish"])
        self.assertTrue(endpoint["uncut_exit_not_claimed_zero"])

    def test_scope_and_correction(self):
        row = q196.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["correction"]["Q195_frontier_skipped_the_immediate_B_source_stage"])
        self.assertFalse(row["scope"]["first_sourced_complex_covariance_independently_rederived"])

    def test_invalid_inputs(self):
        for call in (q196.statuses, q196.endpoint_contract, q196.frontier, q196.review, q196.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q196.rates(0)
        with self.assertRaises(ValueError):
            q196.linear_ledger(0)
        with self.assertRaises(ValueError):
            q196.absorption(-1)
        with self.assertRaises(ValueError):
            q196.quotient_contract(-1)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 445)
        self.assertGreater(record["closure"]["edges"], 1589)
        self.assertGreater(record["closure"]["source_bindings"], 6023)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
