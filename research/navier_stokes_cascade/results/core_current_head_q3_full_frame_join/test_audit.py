from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

import audit
from full_frame import DEPENDENCIES, frame_budget, load_closure, projected_pressure, review


class CurrentHeadQ3FullFrameTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_provider_graph_and_all_hashes(self):
        records, graph, dependencies, sources = load_closure()
        row = review()
        self.assertGreater(len(records), 240)
        self.assertEqual(len(graph), row["graph"]["nodes"])
        self.assertTrue(row["checks"]["provider_graph_is_acyclic"])
        self.assertTrue(all(dependencies.values()))
        self.assertTrue(all(sources.values()))

    def test_all_direct_Q3_provider_statuses(self):
        row = review()
        self.assertEqual(len(DEPENDENCIES), 17)
        self.assertEqual(len(row["provider_status_checks"]), 17)
        self.assertTrue(all(row["provider_status_checks"].values()))

    def test_frame_budget_is_finite_on_closed_path(self):
        row = frame_budget(Fraction(1, 13), Fraction(1, 17), Fraction(1, 19), Fraction(1, 7), 3, 2, Fraction(5, 3), Fraction(7, 4), 13, 17)
        self.assertGreater(row["normal_reciprocal"], 0)
        self.assertGreater(row["tangent_constant"], 0)
        self.assertGreater(row["dstar"], 0)
        self.assertLessEqual(row["picard_next_ratio"], Fraction(1, 2))

    def test_invalid_frame_denominator_or_path_is_rejected(self):
        with self.assertRaises(ValueError):
            frame_budget(0, 1, 1, 1, 1, 0, 1, 1, 1, 1)
        with self.assertRaises(ValueError):
            frame_budget(1, 1, 1, 2, 1, 0, 1, 1, 1, 1)
        with self.assertRaises(ValueError):
            frame_budget(1, 1, 1, 1, 1, 0, 1, 1, 2, 1)

    def test_unprojected_normal_source_pressure_is_retained(self):
        row = projected_pressure((0, 0, 1), (0, 0, 0), (0, 0, 2), (1, 0, 0), 5, 1)
        self.assertEqual(row["source_normal_over_i"], Fraction(2, 5))
        self.assertNotEqual(row["total_without_K_over_i"], 0)

    def test_moving_normal_pressure_is_retained(self):
        row = projected_pressure((2, 1, 1), (1, -1, 0), (0, 0, 0), (1, -2, 0), 7, -2)
        self.assertNotEqual(row["moving_normal_over_i"], 0)

    def test_old_numeric_gate_and_unknown_first_band_are_excluded(self):
        compactness = review()["interface_controls"]["compactness"]
        self.assertFalse(compactness["old_numeric_Cphi_gate_required"])
        self.assertFalse(compactness["unspecified_historical_first_band_certified"])
        self.assertIn("Q2 later tail", compactness["domain"])

    def test_pressure_and_potential_cutoff_commutators_are_recorded(self):
        cutoffs = review()["interface_controls"]["cutoffs"]
        self.assertIn("grad(psi)xA", cutoffs)
        self.assertIn("pi*grad(psi)", cutoffs)

    def test_scope_advances_only_Q3(self):
        scope = review()["scope"]
        self.assertTrue(scope["Q1_Q2_object_is_a_direct_dependency"])
        self.assertTrue(scope["Q3_full_frame_primitive_path_and_phase_domain_bound_on_same_later_tail"])
        self.assertFalse(scope["fixed_projector_compensator_and_weighted_inverse_contracts_closed"])
        self.assertFalse(scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["fixed_projector_compensator_and_weighted_inverse_contracts_closed"] = True
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
