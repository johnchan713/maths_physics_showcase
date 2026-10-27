from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import audit
from provider_reconstruction import (
    EDGES, NODES, compact_fixed_order_bound, compact_remainder_bound, flat_edge_certificate,
    oriented_edge_review, review, slow_response_bound, source_frontier,
    topological_order,
)


class DirectCompactSlowProviderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()

    def test_provider_graph_is_finite_and_acyclic(self):
        self.assertEqual(topological_order(), list(NODES))
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + ((NODES[-1], NODES[0]),))

    def test_flat_weight_absorbs_each_fixed_pole(self):
        for power in (0, 7, 200003):
            row = flat_edge_certificate(power, 1)
            self.assertEqual(row["flat_limit"], 0)
            self.assertTrue(row["fixed_power_absorbed_by_positive_quadratic_exponential"])
        with self.assertRaises(ValueError):
            flat_edge_certificate(7, 0)

    def test_shifted_integrals_have_the_only_valid_edge_orientation(self):
        row = oriented_edge_review()
        self.assertIn("I_-", row["inner_edge"]["active_integral"])
        self.assertIn("I_+", row["outer_edge"]["active_integral"])
        self.assertTrue(row["shift_is_retained"])
        self.assertIn("nonzero total integral", row["wrong_orientation_fails"])

    def test_compact_derivatives_do_not_create_shift_powers(self):
        row = compact_fixed_order_bound(29, 41)
        self.assertEqual(row["input_order"], 29)
        self.assertEqual(row["output_pole"], 41)
        self.assertTrue(row["fixed_coefficient_derivatives_commute_under_shifted_integral"])

    def test_compact_remainder_gain_keeps_zero_mode_condition(self):
        row = compact_remainder_bound(11, 19)
        self.assertEqual(row["source_order"], 52)
        self.assertEqual(row["epsilon_gain"], "19*kappa")
        self.assertTrue(row["zero_haar_weighted_moment_required"])
        with self.assertRaises(ValueError):
            compact_remainder_bound(11, 0)

    def test_slow_response_requests_and_pressure_are_complete(self):
        row = slow_response_bound(12)
        self.assertEqual(row["P_and_Jz_order"], 13)
        self.assertEqual(row["Jtheta_order"], 14)
        self.assertTrue(row["full_axial_velocity_retained"])
        self.assertIn("(partial_Z c_rho)*P", row["moving_pressure_product"])

    def test_operator_bounds_do_not_manufacture_actual_source_bounds(self):
        row = source_frontier()
        self.assertFalse(row["pointwise_band_smoothness_implies_tail_supremum"])
        self.assertFalse(row["operator_estimate_implies_source_estimate"])
        self.assertIn("sup_{bands", row["normalized_source_supremum"])

    def test_complete_residual_scope_is_preserved(self):
        self.assertEqual(set(self.result["residual_ledger"]["groups"]), {
            "time", "viscosity", "pressure", "old_new", "new_old", "new_square",
        })
        self.assertFalse(self.result["scope"]["unconditional_complete_residual_force_comparison_or_blowup_claimed"])

    def test_frozen_record_reproduces_and_tampering_fails(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["remaining_source_frontier"]["proved_from_raw_generated_source_formulas_here"] = True
        with self.assertRaises(ValueError):
            audit.validate(bad)

    def test_direct_dependencies_are_frozen(self):
        record = audit.generate()
        self.assertEqual(len(record["dependency_evidence_sha256"]), 4)
        self.assertTrue(all(record["gates"].values()))

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
