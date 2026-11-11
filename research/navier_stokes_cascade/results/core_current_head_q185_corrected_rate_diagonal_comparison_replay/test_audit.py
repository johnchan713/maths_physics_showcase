import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q185


class Q185Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q185.statuses().values()))

    def test_exact_margin_correction(self):
        row = q185.residual_margin_replay()
        self.assertEqual(row["log_regular_and_cross_margin"], Q(49997, 40000000))
        self.assertEqual(row["pure_regular_and_cross_margin"], Q(49996, 40000000))
        self.assertEqual(row["pure_margin_correction"], -Q(1, 40000000))

    def test_compact_margins_do_not_change(self):
        row = q185.residual_margin_replay()
        self.assertEqual(row["A0_own_margin_both_routes"], Q(1, 40000000))
        self.assertEqual(row["A1_own_margin_both_routes"], Q(1, 40000000))
        self.assertTrue(row["compact_accounts_disjoint"])

    def test_stage_selection_is_exact_and_finite(self):
        for route in ("log", "pure"):
            row = q185.comparison_stage(7, 11, 13, route)
            self.assertTrue(row["threshold_holds"])
            self.assertTrue(row["replay_beats_target_by_margin"])
            self.assertGreaterEqual(row["J0"], 0)

    def test_pure_stage_never_precedes_log_stage(self):
        for request in ((0, 0, 0), (5, 2, 9), (100, 7, 13)):
            log = q185.comparison_stage(*request, route="log")
            pure = q185.comparison_stage(*request, route="pure")
            self.assertGreaterEqual(pure["J0"], log["J0"])

    def test_diagonal_complete_memory_and_object_limit(self):
        row = q185.diagonal_replay()
        self.assertTrue(row["all_six_groups_retain_complete_memory"])
        self.assertEqual(row["localized_term_count"], 12)
        self.assertTrue(row["pure_route_reselects_finite_stages_and_may_change_diagonal_scales"])
        self.assertFalse(row["pure_route_identified_byte_for_byte_with_Q181_field"])

    def test_comparison_and_core_contract(self):
        row = q185.comparison_replay()
        self.assertEqual(row["worst_radius_decay"], 2)
        self.assertEqual(row["cutoff_power"], 8)
        self.assertTrue(row["annular_diagonal_terms_vanish_exactly"])
        self.assertFalse(row["background_terms_vanish_exactly"])
        self.assertFalse(row["unconditional_blowup_claimed"])

    def test_review_scope(self):
        row = q185.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q180_Q182_log_route_replayed_conditionally"])
        self.assertTrue(row["scope"]["Q180_Q182_corrected_pure_route_replayed_conditionally"])
        self.assertFalse(row["scope"]["original_old_rate_exact_pure_power_contract_proved"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            q185.statuses(identity="wrong")
        with self.assertRaises(ValueError):
            q185.comparison_stage(0, route="other")
        for bad in (-1, True, 0.5):
            with self.assertRaises(ValueError):
                q185.comparison_stage(bad)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 431)
        self.assertGreater(record["closure"]["edges"], 1500)
        self.assertGreater(record["closure"]["source_bindings"], 5790)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
