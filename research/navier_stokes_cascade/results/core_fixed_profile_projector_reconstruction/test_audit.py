from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from projectors import (
    EDGES, NODES, arbitrary_family_obstruction, diagnostic_profile,
    dilation_moment, fixed_profile_jet_bound, review, strict_support_inclusion, support_row,
    topological_order,
)


class FixedProfileProjectorJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_dependency_route_is_acyclic(self):
        self.assertEqual(set(topological_order()), set(NODES))

    def test_cycle_is_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + (("fixed_profile_projector_contract", "current_same_tail_full_path"),))

    def test_polynomial_fixture_has_exact_mass_and_A2(self):
        row = diagnostic_profile()
        self.assertEqual(row["mass"], 1)
        self.assertEqual(row["A2"], Q(10, 7))
        self.assertGreater(row["A2"], 0)

    def test_both_projector_moments_are_scale_invariant(self):
        for e in (1, 2):
            for s in (Q(1, 4), Q(1), Q(4)):
                self.assertEqual(dilation_moment(s, e), 1)

    def test_moving_support_uses_actual_scaled_hull(self):
        self.assertEqual(support_row(Q(1, 4)), (Q(1, 2), Q(1)))
        self.assertEqual(support_row(Q(4)), (Q(2), Q(4)))

    def test_actual_support_is_strictly_inside_enclosure(self):
        for scale in (Q(1, 4), Q(1), Q(4)):
            row = strict_support_inclusion(scale, Q(5, 4), Q(7, 4))
            self.assertTrue(row["strict"])
            self.assertNotEqual(row["actual"], row["hull"])
        with self.assertRaises(ValueError):
            strict_support_inclusion(Q(1), Q(1), Q(7, 4))

    def test_normalization_only_family_has_unbounded_Z_jet(self):
        for e in (1, 2):
            one = arbitrary_family_obstruction(e, 1)
            eleven = arbitrary_family_obstruction(e, 11)
            self.assertEqual(one["weighted_derivative_moment"], 0)
            self.assertNotEqual(one["point_derivative"], 0)
            self.assertEqual(eleven["point_derivative"], 11 * one["point_derivative"])

    def test_fixed_profile_majorant_requires_complete_exact_data(self):
        self.assertGreater(fixed_profile_jet_bound(2, [1, 2, 3], [1, 4, 16]), 0)
        for args in ((-1, [], []), (1, [1], [1, 2]), (0, [-1], [1]), (0, [1.0], [1])):
            with self.assertRaises(ValueError):
                fixed_profile_jet_bound(*args)

    def test_same_choice_recompute_is_ordered_after_projector_choice(self):
        order = topological_order()
        self.assertLess(order.index("rho_squared_projectors"), order.index("same_choice_signed_recompute"))
        self.assertLess(order.index("same_choice_signed_recompute"), order.index("fixed_profile_projector_contract"))

    def test_all_controls_pass(self):
        self.assertTrue(self.result["checks"])
        self.assertTrue(all(self.result["checks"].values()))

    def test_dependencies_and_scope(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 9)
        self.assertTrue(all(self.record["gates"].values()))
        flags = self.record["flags"]
        self.assertTrue(flags["normalization_only_projector_jet_inference_rejected"])
        self.assertTrue(flags["current_full_path_tail_used"])
        self.assertTrue(flags["actual_support_strictly_contained_in_moving_enclosure"])
        self.assertTrue(flags["same_choice_source_signed_target_covariance_and_pressure_recomputation_required"])
        for key in (
            "moving_enclosure_identified_with_actual_support",
            "historical_unspecified_projectors_or_targets_identified_with_new_choice",
            "uniformity_over_all_profiles_orders_or_cycles_proved",
            "compact_or_slow_weighted_inverse_domains_closed_here",
            "actual_high_jet_current_fields_or_complete_residual_certified",
            "unconditional_pulse_force_residual_improvement_or_blowup_proved",
        ):
            self.assertFalse(flags[key])

    def test_frozen_record_and_scope_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        changed = deepcopy(frozen)
        changed["flags"]["historical_unspecified_projectors_or_targets_identified_with_new_choice"] = True
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
