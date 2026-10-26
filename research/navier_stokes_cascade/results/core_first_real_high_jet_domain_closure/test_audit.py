from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from high_jet_join import normalized_join, physical_join, primitive_join, review


class FirstRealHighJetDomainClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_all_five_primitive_contracts_are_joined(self):
        row = primitive_join()
        self.assertEqual(len(row), 5)
        self.assertTrue(all(item["independently_joined_on_current_tail"] for item in row.values()))

    def test_normalized_request_indices(self):
        row = normalized_join()
        self.assertEqual(row["Phi_order"], 199994)
        self.assertEqual(row["DzPhi_order"], 199993)
        self.assertEqual(row["Phi_graph_ceiling"], 200018)
        self.assertEqual(row["radial_source_order"], 399995)
        self.assertEqual(row["combined_graph_ceiling"], 400009)

    def test_both_compact_tails_have_strict_kappa_margin(self):
        row = normalized_join()
        self.assertEqual(row["A1_tail"], Q(3) - 5 * row["kappa"])
        self.assertEqual(row["A0_tail"], Q(3) - 5 * row["kappa"])
        self.assertEqual(row["A1_margin"], row["kappa"])
        self.assertEqual(row["A0_margin"], row["kappa"])

    def test_no_compact_output_is_used_as_an_input(self):
        row = normalized_join()
        self.assertTrue(row["all_requests_finite"])
        self.assertFalse(row["desired_compact_outputs_used_as_inputs"])
        self.assertFalse(row["numerical_materialization_required"])

    def test_physical_example_is_finite_and_strict(self):
        row = physical_join()
        self.assertEqual(row["A0_source_order"], 60200007)
        self.assertEqual(row["A1_DzPhi_source_order"], 60000005)
        self.assertEqual(row["combined_graph_ceiling"], 60200021)
        self.assertTrue(row["finite_graph"])
        self.assertGreater(row["A0_margin"], 0)
        self.assertGreater(row["A1_margin"], 0)
        self.assertFalse(row["uniform_in_all_orders_stages_or_scales"])

    def test_invalid_orders_and_parameters_are_rejected(self):
        for value in (-1, True, Q(1, 2)):
            with self.assertRaises(ValueError):
                normalized_join(value)
        for kappa in (0, Q(-1, 100), Q(1, 200000), Q(1, 10)):
            with self.assertRaises(ValueError):
                normalized_join(kappa=kappa)

    def test_controls_dependencies_and_scope(self):
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 9)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertTrue(self.record["flags"]["actual_current_radial_source_fixed_order_inequality_proved_on_certified_tail"])
        self.assertFalse(self.record["flags"]["post_sourced_changed_state_compact_DAG_proved"])
        self.assertFalse(self.record["flags"]["complete_regular_or_physical_residual_improvement_proved"])

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["complete_regular_or_physical_residual_improvement_proved"] = True
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
