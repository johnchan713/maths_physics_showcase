from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from complete_residual import (
    PHYSICAL_MEMORY, RESIDUAL_GROUPS, join, review, shared_source_record,
)


class CurrentHeadQ6ResidualTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_Q5_and_residual_use_identical_source_record(self):
        binding = shared_source_record()
        self.assertTrue(binding["same_digest"])
        self.assertEqual(binding["q5_digest"], binding["residual_digest"])

    def test_default_stage_and_actual_minimum(self):
        row = join()
        self.assertEqual(row["base_stage"], 12816)
        self.assertEqual(row["log_reserved_stage"], 12817)
        self.assertEqual(Q(row["complete_minimum_exponent"]), Q(20030773, 20000000))
        self.assertGreater(Q(row["complete_minimum_exponent"]), 1)

    def test_compact_integration_cost_is_carried_at_same_stage(self):
        row = join(1, 3, (7, 11), 13, 17)
        request = row["stage_request"]
        self.assertEqual(request["stage"], row["log_reserved_stage"])
        for name in ("A0", "A1"):
            self.assertEqual(
                request["raw_sources_with_extra_integration"][name],
                request["raw_sources_before_extra_integration"][name] + 2,
            )

    def test_four_accounts_have_three_disjoint_classes(self):
        accounts = join(1, 3, (7, 11), 13, 17)["reserve_accounts"]
        self.assertEqual(set(accounts), {"A0", "A1", "regular", "cutoff"})
        self.assertEqual(len({entry["owner"] for entry in accounts.values()}), 4)
        self.assertEqual(
            {entry["reserve_class"] for entry in accounts.values()},
            {"compact_source", "regular_residual", "cutoff_tail"},
        )
        self.assertTrue(all(not entry["may_pay_other_account"] for entry in accounts.values()))

    def test_all_six_groups_and_two_crosses_are_retained(self):
        groups = join(1, 3, (7, 11), 13, 17)["six_group_exponents"]
        self.assertEqual(tuple(groups), RESIDUAL_GROUPS)
        self.assertEqual(groups["old_new"], groups["new_old"])
        self.assertIn("new_square", groups)

    def test_common_exponent_is_minimum_not_sum_of_reserves(self):
        row = join(1, 3, (7, 11), 13, 17)
        values = [Q(row["fixed_state_exponent"])]
        values.extend(Q(value) for value in row["six_group_exponents"].values())
        self.assertEqual(Q(row["complete_minimum_exponent"]), min(values))
        self.assertLess(Q(row["complete_minimum_exponent"]), Q(row["fixed_state_exponent"]))

    def test_complete_memory_and_cutoff_terms_are_retained(self):
        self.assertEqual(tuple(self.result["physical_memory"]), PHYSICAL_MEMORY)
        self.assertIn("potential_cutoff_curl_commutators", PHYSICAL_MEMORY)
        self.assertIn("full_pressure_cutoff_derivatives", PHYSICAL_MEMORY)

    def test_fixed_sample_targets_all_close(self):
        self.assertEqual(len(self.result["samples"]), 3)
        self.assertTrue(all(all(row["checks"].values()) for row in self.result["samples"]))

    def test_invalid_orders_are_rejected(self):
        for args in ((-1, 0), (0, -1), (0.0, 0), (True, 0)):
            with self.assertRaises(ValueError):
                join(*args)

    def test_dependency_closure_and_controls(self):
        self.assertEqual(len(self.result["checks"]), 9)
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 6)
        self.assertTrue(all(self.record["gates"].values()))

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["review"]["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"] = True
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
