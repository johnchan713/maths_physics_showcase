from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

import audit
from exact_map import DEPENDENCIES, load_closure, quadratic_root_control, review


class CurrentHeadQ1ExactMapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_provider_graph_and_all_hashes(self):
        records, graph, dependencies, sources = load_closure()
        row = review()
        self.assertGreater(len(records), 238)
        self.assertEqual(len(graph), row["graph"]["nodes"])
        self.assertTrue(row["checks"]["provider_graph_is_acyclic"])
        self.assertTrue(all(dependencies.values()))
        self.assertTrue(all(sources.values()))

    def test_all_direct_Q1_provider_statuses(self):
        row = review()
        self.assertEqual(len(DEPENDENCIES), 14)
        self.assertEqual(len(row["provider_status_checks"]), 10)
        self.assertTrue(all(row["provider_status_checks"].values()))

    def test_annular_and_compact_root_inequalities_are_strict(self):
        annular = quadratic_root_control(610, 50, Fraction(286, 10**19))
        compact = quadratic_root_control(800, 25000, Fraction(1, 2**256))
        for control in (annular, compact):
            self.assertLess(control["invariant"], 1)
            self.assertLess(control["lipschitz"], Fraction(1, 2))
            self.assertGreater(control["radius"], 0)

    def test_pressure_endpoint_obstruction_is_retained(self):
        row = review()["interface_controls"]["pressure"]
        self.assertEqual(row["homogeneous_mode"], "C*(1-eta^2)^(-2A)")
        self.assertFalse(row["nonzero_mode_is_bounded_at_endpoints"])
        self.assertTrue(row["bounded_smooth_solution_forces_C_zero"])

    def test_leading_pressure_does_not_collapse_full_pressure(self):
        text = review()["contracts"]["leading_pressure"]
        for term in ("positive-order", "wave", "mean", "localization"):
            self.assertIn(term, text)

    def test_scope_advances_only_Q1(self):
        scope = review()["scope"]
        self.assertTrue(scope["Q1_provider_chain_bound_on_one_current_head_object"])
        self.assertTrue(scope["Q1_interface_uniqueness_and_pressure_compatibility_rechecked"])
        self.assertFalse(scope["every_provider_theorem_independently_reproved_here"])
        self.assertFalse(scope["Q2_regional_continuum_infima_discharged"])
        self.assertFalse(scope["Q3_full_frame_primitive_domain_discharged"])
        self.assertFalse(scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["Q2_regional_continuum_infima_discharged"] = True
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
