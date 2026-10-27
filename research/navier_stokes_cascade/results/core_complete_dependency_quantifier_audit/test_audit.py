from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import audit
from dependency_graph import ROOTS, load_graph, quantifier_review, review


class CompleteDependencyQuantifierAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_transitive_graph_is_frozen_and_acyclic(self):
        result = review()
        self.assertEqual(result["graph"]["nodes"], 215)
        self.assertEqual(result["graph"]["edges"], 837)
        self.assertTrue(result["checks"]["transitive_evidence_graph_is_acyclic"])

    def test_every_declared_hash_and_source_matches(self):
        _, _, hashes, sources = load_graph()
        self.assertTrue(hashes and all(hashes.values()))
        self.assertTrue(sources and all(sources.values()))

    def test_active_same_object_spine_is_connected(self):
        links = review()["graph"]["spine_links"]
        self.assertTrue(links and all(links.values()))

    def test_fixed_request_quantifier_order(self):
        row = quantifier_review()
        self.assertEqual(row["sample_fixed_comparison_stages"]["m=0,N=0"], 12816)
        self.assertTrue(row["fixed_m_N_choose_one_finite_J_before_q_limit"])
        self.assertTrue(row["comparison_order_is_fixed_T_then_R_then_tau"])

    def test_request_archives_do_not_reselect_the_field(self):
        row = quantifier_review()
        self.assertTrue(row["derivative_and_pole_indices_are_proof_archives_not_field_reselections"])
        self.assertTrue(row["field_and_sequences_chosen_before_requests"])

    def test_provider_frontier_and_scope_remain_explicit(self):
        result = review()
        self.assertFalse(result["quantifiers"]["unconditional_provider_estimates_established_here"])
        self.assertFalse(result["scope"]["unconditional_pulse_residual_force_or_blowup_claimed"])
        self.assertIn("same_fixed_representatives", result["remaining_master_condition"])

    def test_frozen_record_reproduces_and_tampering_fails(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["quantifiers"]["unconditional_provider_estimates_established_here"] = True
        with self.assertRaises(ValueError):
            audit.validate(bad)

    def test_direct_dependencies_are_complete(self):
        self.assertEqual(len(ROOTS), 13)
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 13)
        self.assertTrue(all(self.record["gates"].values()))

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
