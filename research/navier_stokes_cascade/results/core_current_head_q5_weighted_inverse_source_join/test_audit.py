from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

import audit
from weighted_inverse import (
    DEPENDENCIES, fixed_stage_domain, load_closure, log_reserve,
    pure_power_equality, request_graph, review, two_edge_constant,
)


class CurrentHeadQ5WeightedInverseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_provider_graph_and_all_hashes(self):
        records, graph, dependencies, sources = load_closure()
        row = review()
        self.assertGreater(len(records), 242)
        self.assertEqual(len(graph), row["graph"]["nodes"])
        self.assertTrue(row["checks"]["provider_graph_is_acyclic"])
        self.assertTrue(all(dependencies.values()))
        self.assertTrue(all(sources.values()))

    def test_all_direct_Q5_provider_statuses(self):
        row = review()
        self.assertEqual(len(DEPENDENCIES), 10)
        self.assertEqual(len(row["provider_status_checks"]), 10)
        self.assertTrue(all(row["provider_status_checks"].values()))

    def test_three_compact_kernel_branches_are_positive(self):
        row = two_edge_constant(Fraction(4, 49), Fraction(401, 100), Fraction(1, 4900))
        self.assertEqual(row["sum"], row["inner"] + row["outer"] + row["interior"])
        with self.assertRaises(ValueError):
            two_edge_constant(1, 0, 1)

    def test_fast_slow_and_pressure_request_graph(self):
        row = request_graph(0, 399995, 200001)
        self.assertEqual(row["fast_tangential_or_axial"], 4)
        self.assertEqual(row["fast_radial_potential"], 5)
        self.assertEqual(row["prefast_source"], 8)
        self.assertEqual(row["whole_axial_pressure_moment"], 9)
        self.assertEqual(row["A0_with_extra_integration"], 399997)
        self.assertEqual(row["A1_with_extra_integration"], 200003)

    def test_finite_log_degree_is_absorbed_with_half_reserve(self):
        row = log_reserve(Fraction(1, 200), Fraction(1, 100000), 37)
        self.assertGreater(row["retained"], 0)
        self.assertEqual(row["gain"], row["spent_on_logs"] + row["retained"])
        with self.assertRaises(ValueError):
            log_reserve(1, 0, 1)

    def test_old_pure_power_equality_is_preserved_but_reclassified(self):
        row = pure_power_equality(Fraction(1, 100000))
        self.assertTrue(row["equality"])
        self.assertEqual(row["radial_viscosity"], row["nu"])
        self.assertEqual(row["source_margin"] + row["viscosity_margin"], 0)

    def test_fixed_stage_quantifier_order(self):
        row = fixed_stage_domain(7, 13, 19, 37)
        self.assertTrue(row["stage_fixed_first"])
        self.assertFalse(row["stage_uniform_constant_required"])
        self.assertFalse(row["derivative_uniform_constant_required"])

    def test_pressure_and_actual_A1_are_retained(self):
        controls = review()["interface_controls"]
        self.assertIn("gamma=gd-A1*gd", controls["fast_map"])
        self.assertIn("c_rho,Z*P", controls["pressure"])
        self.assertIn("six covariances", controls["pressure"])

    def test_scope_advances_Q5_and_corrects_Q4_pure_power_wording(self):
        scope = review()["scope"]
        self.assertTrue(scope["Q1_Q4_object_is_a_direct_dependency"])
        self.assertTrue(scope["Q5_weighted_inverse_and_log_source_contract_bound_on_same_later_tail"])
        self.assertTrue(scope["finite_stage_log_source_induction_instantiated"])
        self.assertFalse(scope["raw_pure_power_supremum_promoted_as_exact_current_source_domain"])
        self.assertFalse(scope["six_group_complete_residual_replayed"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["six_group_complete_residual_replayed"] = True
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
