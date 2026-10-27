from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
import audit
from current_frontier import DISCHARGE_CHAIN, FRONTIER, load_graph, review

class CurrentHeadConditionalFrontierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.record = audit.generate()
    def test_transitive_graph_counts_and_acyclicity(self):
        result = review(); self.assertEqual(result["graph"]["nodes"], 224); self.assertEqual(result["graph"]["edges"], 893); self.assertTrue(result["checks"]["current_head_graph_is_acyclic"])
    def test_every_declared_binding_matches(self):
        _, _, dependencies, sources = load_graph(); self.assertEqual(len(dependencies), 893); self.assertEqual(len(sources), 3675); self.assertTrue(all(dependencies.values())); self.assertTrue(all(sources.values()))
    def test_discharge_chain_is_connected(self):
        links = review()["graph"]["chain_links"]; self.assertEqual(len(links), len(DISCHARGE_CHAIN) - 1); self.assertTrue(all(links.values()))
    def test_five_contract_frontier_is_live(self):
        result = review(); self.assertEqual(len(result["live_frontier"]), len(FRONTIER)); self.assertTrue(all(result["graph"]["frontier_links"].values()))
    def test_old_provider_condition_is_not_reported_as_live_frontier(self):
        row = review(); self.assertIn("actual_source_induction", row["discharged_relative_to_frontier"]); self.assertIn("complete_residual", row["discharged_relative_to_frontier"]); self.assertNotIn("remaining_master_condition", row)
    def test_same_object_and_limit_guards(self): self.assertTrue(all(review()["quantifier_guard"].values()))
    def test_scope_remains_conditional(self):
        scope = review()["scope"]; self.assertTrue(scope["downstream_chain_complete_relative_to_five_contracts"]); self.assertFalse(scope["five_contracts_independently_reproved_in_this_checkpoint"]); self.assertFalse(scope["hashes_or_tests_treated_as_continuum_proof"]); self.assertFalse(scope["unconditional_pulse_force_exclusion_or_blowup_claimed"])
    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text()); audit.validate(frozen); bad = deepcopy(frozen); bad["review"]["scope"]["unconditional_pulse_force_exclusion_or_blowup_claimed"] = True
        with self.assertRaises(ValueError): audit.validate(bad)
    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"; path.write_text("keep")
            with self.assertRaises(SystemExit): audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")
if __name__ == "__main__": unittest.main()
