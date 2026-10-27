from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from compactness import (
    finite_dag_majorant, implicit_floor, primitive_ledger,
    reciprocal_jet_bound, review, stage_row,
)


class ClosedChartCoefficientJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_implicit_floor_includes_both_angular_poles(self):
        self.assertEqual(implicit_floor(Q(1, 100)), Q(49, 50))

    def test_reciprocal_jet_recurrence_is_finite_at_every_fixed_order(self):
        rows = reciprocal_jet_bound(Q(49, 50), 4, 12)
        self.assertEqual(rows[0], Q(50, 49))
        self.assertTrue(all(x > 0 for x in rows))

    def test_every_primitive_contract_uses_a_closed_domain(self):
        ledger = primitive_ledger()
        self.assertEqual(len(ledger), 7)
        self.assertTrue(all(row["closed"] for row in ledger.values()))

    def test_fixed_stage_graphs_are_finite_and_support_preserving(self):
        for args in ((0, 0, 0), (1, 2, 1), (2, 3, 5)):
            row = stage_row(*args)
            self.assertTrue(row["finite_and_acyclic"])
            self.assertTrue(row["closed_supports_preserved"])
            self.assertGreaterEqual(row["physical_primitive_ceiling"], row["state_primitive_ceiling"])

    def test_finite_composition_proxy_never_uses_an_infinite_graph(self):
        values = [finite_dag_majorant(2, Q(1, 10), n) for n in range(5)]
        self.assertTrue(all(x > 0 for x in values))
        self.assertEqual(values, sorted(values))

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError): implicit_floor(0)
        with self.assertRaises(ValueError): reciprocal_jet_bound(0, 1, 1)
        with self.assertRaises(ValueError): finite_dag_majorant(1, 0, 1)
        with self.assertRaises(ValueError): stage_row(-1)

    def test_controls_dependencies_and_conditional_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 7)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))
        self.assertFalse(result["scope"]["actual_exact_core_background_and_primitive_margins_independently_proved_here"])

    def test_frozen_record_reproduces_and_scope_tampering_fails(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["unconditional_complete_residual_compact_force_or_blowup"] = True
        with self.assertRaises(ValueError): audit.validate(bad)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit): audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
