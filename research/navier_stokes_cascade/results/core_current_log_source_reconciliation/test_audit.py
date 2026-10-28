from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from reconcile import (
    EDGES, MEMORY, NODES, current_leaf_domain, fixed_request_bridge,
    review, topological_order,
)


class CurrentLogSourceReconciliationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_current_chain_is_acyclic_and_ordered(self):
        self.assertEqual(topological_order(), list(NODES))

    def test_cycle_and_unknown_node_are_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + ((NODES[-1], NODES[0]),))
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + (("missing", NODES[0]),))

    def test_fixed_current_leaves_have_log_degree_zero(self):
        row = current_leaf_domain()
        self.assertTrue(all(value == 0 for key, value in row.items() if key != "interpretation"))

    def test_invalid_fixed_requests_are_rejected(self):
        for args in ((-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1), (0.0, 0, 0, 0)):
            with self.assertRaises(ValueError):
                fixed_request_bridge(*args)

    def test_extra_integration_cost_is_two_derivatives(self):
        row = fixed_request_bridge(1, 3, 2, 7)
        for name in ("A0", "A1"):
            self.assertEqual(
                row["raw_sources_with_extra_integration"][name],
                row["raw_sources_before_extra_integration"][name] + 2,
            )
        self.assertEqual(row["extra_source_derivatives_per_compact_map"], 2)

    def test_fixed_log_absorption_keeps_positive_margin(self):
        row = fixed_request_bridge(2, 7, 5, 37)
        self.assertGreater(Q(row["absorption_exponent_bound"]), 0)
        self.assertTrue(all(Q(value) > 0 for value in row["retained_strict_margin"].values()))

    def test_large_finite_request_is_accepted_without_uniformity(self):
        row = fixed_request_bridge(5, 11, 13, 113)
        self.assertTrue(row["finite_formula_graph"])
        self.assertTrue(row["current_inverse_accepts_every_fixed_order"])
        self.assertGreater(row["raw_sources_with_extra_integration"]["A0"], 10**9)

    def test_complete_memory_is_retained(self):
        self.assertEqual(tuple(self.result["memory"]), MEMORY)

    def test_all_controls_and_dependency_sources(self):
        self.assertEqual(len(self.result["checks"]), 12)
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 9)
        self.assertTrue(all(self.record["gates"].values()))

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["complete_residual_or_downstream_replays_reproved_here"] = True
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_frontier_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["pure_power_zero_margin_reintroduced_as_current_frontier"] = True
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
