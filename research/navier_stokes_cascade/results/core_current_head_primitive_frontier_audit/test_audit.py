from copy import deepcopy
import json
import math
from pathlib import Path
import tempfile
import unittest

import audit
from current_frontier import (
    BOUNDARY_WITNESSES, CURRENT_REPLAYS, OLD_FRONTIER, RECONSTRUCTIONS,
    load_graph, review, selector_ceiling,
)


class CurrentHeadPrimitiveFrontierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_transitive_graph_is_larger_than_historical_audit_and_acyclic(self):
        row = review()
        self.assertGreater(row["graph"]["nodes"], 224)
        self.assertGreater(row["graph"]["edges"], 893)
        self.assertTrue(row["checks"]["current_head_graph_is_acyclic"])

    def test_every_declared_binding_matches(self):
        _, _, dependencies, sources = load_graph()
        self.assertEqual(len(dependencies), review()["graph"]["dependency_hash_checks"])
        self.assertEqual(len(sources), review()["graph"]["source_hash_checks"])
        self.assertTrue(all(dependencies.values()))
        self.assertTrue(all(sources.values()))

    def test_six_reconstructions_and_six_current_replays_are_connected(self):
        row = review()
        self.assertEqual(len(row["graph"]["reconstruction_links"]), len(RECONSTRUCTIONS))
        self.assertTrue(all(row["graph"]["reconstruction_links"].values()))
        self.assertEqual(len(row["graph"]["replay_links"]), len(CURRENT_REPLAYS) - 1)
        self.assertTrue(all(row["graph"]["replay_links"].values()))

    def test_old_five_contract_frontier_is_historical(self):
        row = review()
        self.assertFalse(row["correction"]["historical_five_contract_frontier_is_live"])
        self.assertEqual(len(row["graph"]["historical_frontier_links"]), len(OLD_FRONTIER))
        self.assertTrue(all(row["graph"]["historical_frontier_links"].values()))

    def test_lower_frontier_has_three_reached_witnesses(self):
        row = review()
        self.assertEqual(len(row["live_frontier"]), 3)
        self.assertEqual(len(row["graph"]["boundary_links"]), len(BOUNDARY_WITNESSES))
        self.assertTrue(all(row["graph"]["boundary_links"].values()))

    def test_selector_ceiling_requires_actual_positive_continuum_rows(self):
        self.assertEqual(selector_ceiling(0.5, [8.0, 1.0], [1.0, 0.5]), 5)
        for bad in (0.0, -1.0, math.inf, math.nan):
            with self.assertRaises(ValueError):
                selector_ceiling(0.5, [8.0], [bad])
        with self.assertRaises(ValueError):
            selector_ceiling(0.5, [math.inf], [1.0])

    def test_same_object_and_limit_guards(self):
        self.assertTrue(all(review()["quantifier_guard"].values()))

    def test_scope_remains_conditional(self):
        scope = review()["scope"]
        self.assertTrue(scope["former_five_contracts_reconstructed_on_current_objects"])
        self.assertTrue(scope["downstream_replay_complete_relative_to_lower_frontier"])
        self.assertFalse(scope["Q1_Q2_Q3_independently_reproved_here"])
        self.assertFalse(scope["hashes_or_tests_treated_as_continuum_proof"])
        self.assertFalse(scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"] = True
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
