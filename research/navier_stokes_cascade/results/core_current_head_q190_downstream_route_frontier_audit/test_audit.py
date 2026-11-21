import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q190


class Q190Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q190.statuses().values()))

    def test_exact_comparison_stages(self):
        for route, margin in (("log", Q(49997, 40000000)), ("pure", Q(49996, 40000000))):
            for request in ((0, 0, 0), (7, 11, 13), (50, 400, 900)):
                row = q190.comparison_stage(*request, route=route)
                self.assertTrue(row["threshold_holds"])
                self.assertTrue(row["replay_holds_with_margin"])
                self.assertEqual(row["retained_margin"], margin)

    def test_quantifier_order(self):
        row = q190.quantifier_order()
        self.assertTrue(row["all_strict"])
        self.assertFalse(row["comparison_stage_depends_on_q"])
        self.assertFalse(row["stage_uniform_constant_required"])
        self.assertEqual(row["A0_A1_each_own_margin"], Q(1, 40000000))

    def test_downstream_chain(self):
        self.assertTrue(all(q190.downstream_contract().values()))

    def test_frontier_correction_and_scope(self):
        row = q190.review()
        self.assertTrue(all(row["checks"].values()))
        frontier = row["correction_and_frontier"]
        self.assertTrue(frontier["Q185_named_exact_core_and_finite_maps_as_next"])
        self.assertTrue(frontier["Q164_Q165_and_Q189_now_discharge_that_provider_level_frontier"])
        self.assertFalse(frontier["missing_scalar_selector_or_diagonal_inequality_through_Q185"])
        self.assertFalse(frontier["missing_named_analytical_provider_through_Q185"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        for call in (
            q190.statuses, q190.quantifier_order, q190.downstream_contract,
            q190.correction_and_frontier, q190.review, q190.record,
        ):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        for args in ((True, 0, 0), (-1, 0, 0), (0, -1, 0), (0, 0, -1)):
            with self.assertRaises(ValueError):
                q190.comparison_stage(*args)
        with self.assertRaises(ValueError):
            q190.comparison_stage(0, route="other")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 440)
        self.assertGreater(record["closure"]["edges"], 1550)
        self.assertGreater(record["closure"]["source_bindings"], 5970)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
