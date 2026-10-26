from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from frontier import EDGES, NODES, PRIMITIVES, finite_conjunction_threshold, pressure_shape_bound, topological_order


class CorrectedFoundationFrontierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_graph_is_acyclic_and_complete(self):
        order = topological_order()
        self.assertEqual(set(order), set(NODES))
        self.assertEqual(len(order), len(NODES))

    def test_cycle_is_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + (("corrected_comparison", "fixed_point"),))

    def test_finite_tail_is_maximum(self):
        self.assertEqual(finite_conjunction_threshold((17, 4, 29, 11)), 29)
        with self.assertRaises(ValueError):
            finite_conjunction_threshold(())

    def test_pressure_shape_keeps_four_shifts(self):
        b0 = pressure_shape_bound(0, F(1, 2), 2, 3, 5)
        b1 = pressure_shape_bound(1, F(1, 2), 2, 3, 5)
        self.assertEqual(b1, 5 * b0)

    def test_pressure_shape_rejects_bad_data(self):
        with self.assertRaises(ValueError):
            pressure_shape_bound(-1, 1, 1, 1, 1)
        with self.assertRaises(ValueError):
            pressure_shape_bound(0, 1, -1, 1, 1)

    def test_stale_frontier_is_corrected_without_scope_promotion(self):
        flags = self.record["flags"]
        self.assertTrue(flags["leading_core_pressure_is_identified_with_same_completed_backward_pressure"])
        self.assertFalse(flags["corrected_head_next_obligation_text_is_current_frontier"])
        self.assertFalse(flags["five_primitive_accepted_tail_contracts_independently_certified_here"])
        self.assertFalse(flags["unconditional_pulse_residual_force_or_blowup_proved"])

    def test_primitive_frontier_is_explicit(self):
        self.assertEqual(len(PRIMITIVES), 5)
        self.assertEqual(self.record["review"]["primitive_contracts"], list(PRIMITIVES))

    def test_dependencies_and_sources_are_frozen(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 20)
        self.assertTrue(self.record["inherited_source_checks"])
        self.assertTrue(all(self.record["inherited_source_checks"].values()))

    def test_frozen_record_and_scope_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["flags"]["unconditional_pulse_residual_force_or_blowup_proved"] = True
        with self.assertRaises(ValueError):
            audit.validate(bad)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
