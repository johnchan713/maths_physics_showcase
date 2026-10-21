#!/usr/bin/env python3
"""Tests for the foundation-bundle route checkpoint."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from route_checks import EDGES, NODES, chart_domain, combined_tail, review, topological_order


class FoundationBundleRouteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_route_is_acyclic_and_complete(self):
        order = topological_order()
        self.assertEqual(set(order), set(NODES))
        position = {node: index for index, node in enumerate(order)}
        self.assertTrue(all(position[a] < position[b] for a, b in EDGES))

    def test_cycle_is_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(NODES, EDGES + (("complete_residual_force", "core_contraction"),))

    def test_unknown_node_and_duplicate_are_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(NODES, EDGES + (("missing", "core_growth"),))
        with self.assertRaises(ValueError):
            topological_order(NODES + (NODES[0],), EDGES)

    def test_finite_gate_conjunction_is_maximum(self):
        self.assertEqual(combined_tail(16, 31, 17, 29), 31)
        self.assertEqual(combined_tail(10**200, 17), 10**200)
        with self.assertRaises(ValueError):
            combined_tail()
        with self.assertRaises(TypeError):
            combined_tail(1, F(2, 1))

    def test_certified_domain_margins(self):
        domain = chart_domain()
        self.assertEqual(domain["unlocalized_whole_band"], F(1, 4))
        self.assertEqual(domain["localized_q"], F(3, 128))
        self.assertEqual(domain["localized_whole_chart"], F(3, 8))
        self.assertTrue(all(value < 1 for value in domain.values()))

    def test_dependency_routes_and_bindings(self):
        self.assertTrue(all(self.record["reviews"]["dependency_routing"]["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 22)
        self.assertTrue(review()["all_edges_forward"])

    def test_scope_does_not_turn_routing_into_proof(self):
        for key in (
            "historical_first_band_field_or_force_identified",
            "provider_continuum_theorems_independently_reproved_here",
            "frozen_boolean_evidence_treated_as_continuum_proof",
            "formal_or_external_independent_certificate",
            "unconditional_pulse_residual_force_or_blowup_proved",
        ):
            self.assertIs(self.record["flags"][key], False)

    def test_false_scope_promotion_is_rejected(self):
        for key, value in audit.FLAGS.items():
            if not value:
                changed = deepcopy(self.record)
                changed["flags"][key] = True
                with self.assertRaises(ValueError):
                    audit.validate(changed)

    def test_frontier_cannot_be_erased(self):
        changed = deepcopy(self.record)
        changed["remaining_obligations"].pop()
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_frozen_record_and_fresh_replay(self):
        audit.validate(json.loads((audit.HERE / "evidence.json").read_text()))
        with tempfile.TemporaryDirectory(dir=audit.PROJECT) as directory:
            output = Path(directory) / "fresh.json"
            audit.main([
                "--verify-record", str(audit.HERE / "evidence.json"),
                "--output", str(output),
            ])
            audit.validate(json.loads(output.read_text()))
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(output)])


if __name__ == "__main__":
    unittest.main()
