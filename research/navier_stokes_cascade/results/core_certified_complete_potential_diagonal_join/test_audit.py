from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from diagonal_join import (
    certified_input_review, closure_statement, fixed_comparison_review,
    obstruction_correction, review, scale_absorption_review,
)


class CertifiedCompletePotentialDiagonalJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_certified_stage_inputs_are_used(self):
        row = certified_input_review()
        self.assertEqual(row["state_zero_primitive_ceiling"], 900047)
        self.assertEqual(row["physical_zero_primitive_ceiling"], 60900051)
        self.assertEqual(row["stage_two_distinct_poles"], [Q(371400672)])
        self.assertTrue(row["all_fixed_graphs_finite_and_acyclic"])

    def test_arbitrary_finite_constants_are_absorbed(self):
        row = scale_absorption_review()
        self.assertTrue(row["all_extreme_margins_strict"])
        self.assertTrue(all(row["inherited_scale_controls"]["checks"].values()))
        self.assertFalse(row["stage_uniform_constant_bound_used"])
        self.assertFalse(row["raw_series_convergence_used"])

    def test_previous_obstruction_is_scoped_to_naive_substitution(self):
        row = obstruction_correction()
        self.assertFalse(row["naive_uncut_substitution_proved"])
        self.assertFalse(row["naive_stage_growth_bound_available"])
        self.assertFalse(row["naive_inequality_required_by_locally_finite_cutoff_diagonal"])
        self.assertTrue(row["arbitrary_finite_stage_constants_absorbed"])

    def test_fixed_comparison_closes_complete_residual(self):
        row = fixed_comparison_review()
        self.assertEqual(row["comparison_stage"], 12816)
        self.assertTrue(all(row["checks"].values()))
        self.assertEqual(set(row["complete_residual_groups"]), {
            "time", "viscosity", "pressure", "old_new", "new_old", "new_square",
        })
        self.assertTrue(row["comparison_stage_fixed_before_q_limit"])
        self.assertFalse(row["uncontrolled_j_of_q_substitution_used"])

    def test_multiple_fixed_comparisons(self):
        for h in (Q(1, 200), Q(9, 1000), Q(1, 10000)):
            for m, target in ((0, 0), (1, 3), (5, 11), (20, 2)):
                self.assertTrue(all(fixed_comparison_review(m, target, h)["checks"].values()))

    def test_scope_is_conditional_and_global_work_remains(self):
        row = closure_statement()
        self.assertTrue(row["one_locally_finite_divergence_free_field_follows_conditionally"])
        self.assertTrue(row["all_fixed_jets_of_complete_local_residual_are_flat_conditionally"])
        self.assertFalse(row["stage_uniform_constants_or_naive_j_of_q_bound_required"])
        self.assertFalse(row["actual_infinite_coefficients_or_cutoff_scales_evaluated"])
        self.assertFalse(row["closed_time_localization_global_force_energy_comparison_or_blowup_proved"])

    def test_controls_dependencies_and_scope(self):
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertTrue(self.record["flags"]["previous_stage_growth_obstruction_corrected"])
        self.assertFalse(self.record["flags"]["global_force_energy_comparison_or_blowup_proved"])

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["global_force_energy_comparison_or_blowup_proved"] = True
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
