import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q216


class Q216Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q216.statuses().values()))

    def test_borel_split_without_materialization(self):
        row = q216.borel_split(q216.R_STAR)
        self.assertEqual(row["finite_head_count"], q216.R_STAR - 1)
        self.assertEqual(row["geometric_tail_start"], q216.R_STAR)
        self.assertTrue(row["tail_selector_contains_requested_jet"])
        self.assertFalse(row["requires_order_uniform_constant"])
        self.assertEqual(q216.borel_split(0)["finite_head_count"], 0)
        with self.assertRaises(ValueError):
            q216.borel_split(-1)

    def test_source_quantifiers(self):
        row = q216.source_quantifiers()
        self.assertEqual(row["safe_request"], 2_626_513_221_816)
        self.assertTrue(row["background_finite_head_plus_geometric_tail"])
        self.assertTrue(row["infinite_band_frame_uses_analytic_tail_not_compactness"])
        self.assertTrue(row["compact_kernel_accepts_every_finite_pole"])
        self.assertTrue(row["twelve_q_is_terminal_requested_pole"])
        self.assertFalse(row["diagnostic_fixture_used_as_continuum_proof"])
        self.assertFalse(row["circular_dependency_found"])

    def test_selection_order(self):
        row = q216.selection_order()
        self.assertEqual(len(row), 8)
        self.assertEqual(len(row), len(set(row)))
        self.assertTrue(row.index("exact core and selected roots") < row.index("residual constant then diagonal scale"))

    def test_route_reserves(self):
        row = q216.reserves()
        self.assertEqual(row["regular_log_reserve"], Q(49997, 40000000))
        self.assertEqual(row["compact_log_reserve"], Q(1, 40000000))
        self.assertEqual(row["pure_power_margin"], 0)
        self.assertTrue(row["pure_power_route_rejected"])

    def test_review(self):
        row = q216.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q215_provider_quantifiers_supported_by_source_recursions"])
        self.assertFalse(row["scope"]["R_star_numerically_materialized"])
        self.assertFalse(row["scope"]["unconditional_force_exclusion_or_blowup"])

    def test_identity_guards(self):
        for call in (q216.statuses, q216.source_quantifiers, q216.selection_order,
                     q216.reserves, q216.review, q216.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q216.borel_split(4, identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 466)
        self.assertGreaterEqual(record["closure"]["edges"], 1740)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6250)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
