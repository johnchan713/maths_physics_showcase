from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

import audit
from fixed_profiles import (
    DEPENDENCIES, fixed_jet_majorant, load_closure, moving_pressure_memory,
    projector_contract, review,
)


class CurrentHeadQ4FixedProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_provider_graph_and_all_hashes(self):
        records, graph, dependencies, sources = load_closure()
        row = review()
        self.assertGreater(len(records), 241)
        self.assertEqual(len(graph), row["graph"]["nodes"])
        self.assertTrue(row["checks"]["provider_graph_is_acyclic"])
        self.assertTrue(all(dependencies.values()))
        self.assertTrue(all(sources.values()))

    def test_all_direct_Q4_provider_statuses(self):
        row = review()
        self.assertEqual(len(DEPENDENCIES), 10)
        self.assertEqual(len(row["provider_status_checks"]), 10)
        self.assertTrue(all(row["provider_status_checks"].values()))

    def test_fixed_profile_majorant_requires_every_fixed_jet(self):
        self.assertGreater(fixed_jet_majorant(2, [1, 2, 3], [1, 4, 16]), 0)
        for args in ((-1, [], []), (1, [1], [1, 2]), (0, [-1], [1])):
            with self.assertRaises(ValueError):
                fixed_jet_majorant(*args)

    def test_both_exact_projector_moments_and_supports(self):
        for exponent in (1, 2):
            for scale in (Fraction(1, 4), Fraction(1), Fraction(4)):
                row = projector_contract(1, 2, Fraction(5, 4), Fraction(7, 4), scale, Fraction(10, 7), exponent)
                self.assertEqual(row["moment"], 1)
                self.assertTrue(row["strict_support"])
                self.assertNotEqual(row["actual_support"], row["hull"])

    def test_invalid_normalization_or_support_is_rejected(self):
        with self.assertRaises(ValueError):
            projector_contract(1, 2, 1, Fraction(7, 4), 1, Fraction(10, 7), 1)
        with self.assertRaises(ValueError):
            projector_contract(1, 2, Fraction(5, 4), Fraction(7, 4), 1, Fraction(1, 2), 1)

    def test_whole_moving_pressure_memory_is_retained(self):
        row = moving_pressure_memory(Fraction(1, 5), Fraction(7, 3), Fraction(11, 4), Fraction(13, 6), Fraction(17, 9))
        self.assertNotEqual(row["crho_Z_times_P"], 0)
        self.assertEqual(row["total_without_epsilon"], row["Jz_Z"] + row["crho_Z_times_P"] + row["crho_times_P_Z"])

    def test_same_choice_recompute_and_all_six_covariances_recorded(self):
        row = review()["interface_controls"]
        self.assertIn("six covariances", row["same_choice"])
        self.assertIn("c_rho,Z*P", row["commutator"])

    def test_scope_advances_only_Q4(self):
        scope = review()["scope"]
        self.assertTrue(scope["Q1_Q3_field_is_a_direct_dependency"])
        self.assertTrue(scope["Q4_fixed_profile_compensator_projectors_bound_on_same_later_tail"])
        self.assertTrue(scope["same_choice_projector_dependent_objects_recomputed"])
        self.assertFalse(scope["compact_slow_weighted_inverse_contract_closed"])
        self.assertFalse(scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["historical_unspecified_projectors_identified_with_new_choice"] = True
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
