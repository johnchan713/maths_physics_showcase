from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from residual import (
    MEMORY, RESIDUAL_GROUPS, comparison, current_stage_request,
    reserve_accounts, review,
)


class CurrentLogReservedResidualTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_current_stage_request_carries_two_derivative_cost(self):
        row = current_stage_request(17, 3, 5)
        for name in ("A0", "A1"):
            self.assertEqual(
                row["raw_sources_with_extra_integration"][name],
                row["raw_sources_before_extra_integration"][name] + 2,
            )
        self.assertTrue(row["finite_formula_graph"])

    def test_invalid_stage_requests_are_rejected(self):
        for args in ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0.0, 0, 0)):
            with self.assertRaises(ValueError):
                current_stage_request(*args)

    def test_four_named_accounts_use_three_disjoint_classes(self):
        row = comparison(1, 3, (7, 11), 13, 17)
        accounts = row["reserve_accounts"]
        self.assertEqual(set(accounts), {"A0", "A1", "regular", "cutoff"})
        self.assertEqual(len({entry["owner"] for entry in accounts.values()}), 4)
        self.assertEqual(
            {entry["reserve_class"] for entry in accounts.values()},
            {"compact_source", "regular_residual", "cutoff_tail"},
        )
        self.assertTrue(all(not entry["may_pay_other_account"] for entry in accounts.values()))

    def test_account_inputs_are_validated(self):
        base = comparison()
        with self.assertRaises(ValueError):
            reserve_accounts(base, (0,), 0, 0)
        for args in (((-1, 0), 0, 0), ((0, 0), -1, 0), ((0, 0), 0, -1)):
            with self.assertRaises(ValueError):
                reserve_accounts(base, *args)

    def test_default_current_residual_stage_and_target(self):
        row = comparison()
        self.assertEqual(row["base_comparison_stage"], 12816)
        self.assertEqual(row["log_reserved_comparison_stage"], 12817)
        self.assertTrue(all(row["checks"].values()))
        self.assertGreaterEqual(Q(row["complete_minimum_exponent"]), 1)

    def test_large_fixed_stage_request_remains_finite(self):
        row = comparison(5, 11, (37, 41), 43, 47)
        request = row["current_stage_request"]
        self.assertTrue(request["finite_formula_graph"])
        self.assertGreater(request["raw_sources_with_extra_integration"]["A0"], 10**8)

    def test_all_six_groups_and_both_crosses_are_retained(self):
        row = comparison(1, 3, (7, 11), 13, 17)
        groups = row["residual_group_exponents"]
        self.assertEqual(tuple(groups), RESIDUAL_GROUPS)
        self.assertEqual(groups["old_new"], groups["new_old"])
        self.assertIn("Delta e", row["residual_identity"])
        self.assertIn("grad pi", row["residual_identity"])

    def test_complete_memory_is_retained(self):
        self.assertEqual(tuple(self.result["memory"]), MEMORY)

    def test_all_controls_and_dependency_sources(self):
        self.assertEqual(len(self.result["checks"]), 8)
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)
        self.assertTrue(all(self.record["gates"].values()))

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["corner_endpoint_force_energy_or_comparison_replayed_here"] = True
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_claim_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["unconditional_pulse_residual_force_exclusion_or_blowup_proved"] = True
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_frozen_record_and_no_overwrite(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
