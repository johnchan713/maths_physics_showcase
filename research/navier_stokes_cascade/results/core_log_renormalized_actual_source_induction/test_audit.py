from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from source_induction import (
    MEMORY, corrected_source_frontier, fixed_stage_review,
    logarithmic_absorption, no_log_counterexample, one_cycle_review,
    review, seed_review,
)


class LogRenormalizedActualSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()

    def test_pure_power_frontier_is_too_strong(self):
        self.assertEqual(no_log_counterexample(9)["pure_power_ratio"], 10)
        self.assertEqual(no_log_counterexample(1000)["pure_power_ratio"], 1001)

    def test_log_renormalized_counterexample_is_uniform(self):
        for ell in (0, 1, 10, 1000):
            self.assertEqual(no_log_counterexample(ell)["degree_one_renormalized_ratio"], 1)

    def test_positive_slack_absorbs_fixed_log_degree(self):
        row = logarithmic_absorption(11, Q(1, 200000))
        self.assertEqual(Q(row["exponent_bound"]), Q(121 * 200000, 4))
        with self.assertRaises(ValueError):
            logarithmic_absorption(11, 0)

    def test_seed_uses_actual_source_formula(self):
        row = seed_review()
        self.assertEqual(row["source_formula"], "linear + complete_wave_wave + complete_mean_wave")
        self.assertTrue(row["every_prefast_group_meets_mu"])
        self.assertTrue(row["full_pressure_recomputed_after_fast"])

    def test_one_cycle_keeps_actual_sourced_and_generic_rows(self):
        row = one_cycle_review()
        self.assertTrue(row["no_desired_source_callback_used"])
        self.assertEqual(row["memory"], list(MEMORY))
        self.assertEqual(row["generic"]["whole_regular"], row["generic"]["A1_source"])

    def test_fixed_stage_graph_and_archives_are_finite(self):
        for stage in (0, 1, 2, 7):
            row = fixed_stage_review(stage)
            self.assertTrue(row["finite_formula_graph"])
            self.assertGreater(row["first_real_interface_ceiling"], 0)
            self.assertTrue(row["distinct_common_poles"])

    def test_one_extra_integration_pays_for_logs(self):
        row = fixed_stage_review(1, 2, 3)
        self.assertEqual(Q(row["spend_on_logs"]), Q(1, 40000000))
        self.assertTrue(all(Q(v) > 0 for v in row["retained_strict_margin"].values()))

    def test_corrected_frontier_retains_log_degree_and_weight(self):
        row = corrected_source_frontier()
        self.assertIn("(1+|log epsilon_l|)^(-P_(j,R))", row["proved_domain"])
        self.assertIn("zeta^(-1) delta^p", row["proved_domain"])
        self.assertTrue(row["finite_log_degree_for_each_fixed_request"])

    def test_full_memory_and_six_residual_groups_remain(self):
        self.assertEqual(self.result["corrected_source_frontier"]["memory"], list(MEMORY))
        self.assertEqual(set(self.result["residual_groups"]), {
            "time", "viscosity", "pressure", "old_new", "new_old", "new_square",
        })
        self.assertFalse(self.result["scope"]["six_group_complete_residual_replayed_with_log_reserve"])

    def test_frozen_record_reproduces_and_tampering_fails(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["pure_power_source_supremum_without_spending_slack_proved"] = True
        with self.assertRaises(ValueError):
            audit.validate(bad)

    def test_six_dependencies_are_frozen(self):
        record = audit.generate()
        self.assertEqual(len(record["dependency_evidence_sha256"]), 6)
        self.assertTrue(all(record["gates"].values()))

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
