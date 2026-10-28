from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from corner_endpoint import H, KAPPA, corner_target, endpoint_request, review, selected_chart


class CurrentHeadQ7CornerEndpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_selected_chart_uses_exact_Q6_parameters(self):
        chart = selected_chart()
        self.assertEqual(chart["h"], str(H))
        self.assertEqual(chart["kappa"], str(KAPPA))
        self.assertTrue(chart["same_Q6_parameter_choice"])
        self.assertFalse(chart["alternate_h_diagnostic_promoted"])
        self.assertEqual(chart["terminal_identity"], 1)
        self.assertEqual(chart["jacobian_floor"], Q(99, 100))

    def test_default_Q6_stage_reaches_corner(self):
        row = corner_target()
        self.assertEqual(row["comparison_stage"], 12817)
        self.assertEqual(row["stage_request"]["stage"], 12817)
        self.assertGreaterEqual(Q(row["closed_chart_exponent"]), 1)
        self.assertEqual(row["cartesian_factor"], 4)

    def test_extra_compact_derivatives_are_retained(self):
        request = corner_target(1, 3, (7, 11), 13, 17)["stage_request"]
        for name in ("A0", "A1"):
            self.assertEqual(
                request["raw_sources_with_extra_integration"][name],
                request["raw_sources_before_extra_integration"][name] + 2,
            )

    def test_arbitrary_fixed_corner_powers_close(self):
        for jet, power in ((0, 0), (3, 5), (5, 11)):
            row = corner_target(jet, power)
            self.assertTrue(all(row["checks"].values()))
            self.assertGreaterEqual(Q(row["retained_exponent_margin"]), 0)
            self.assertEqual(row["cartesian_factor"], 4 ** (power + 1))

    def test_endpoint_derivative_ledger_is_exact(self):
        for a, j, power in ((0, 0, 0), (1, 3, 5), (5, 4, 11), (3, 7, 2)):
            row = endpoint_request(a, j, power)
            total = a + j
            self.assertEqual(row["punctured_field_Cauchy_orders"], {
                "potential": total + 4,
                "direct_swirl": total + 3,
                "pressure": total + 2,
            })
            self.assertEqual(row["residual_Cauchy_jet"], total + 1)
            self.assertTrue(all(row["checks"].values()))

    def test_higher_request_changes_stage_not_blocking_inequality(self):
        low = endpoint_request(0, 0, 0)
        high = endpoint_request(5, 5, 10)
        self.assertGreater(high["comparison_stage"], low["comparison_stage"])
        self.assertEqual(low["endpoint_blocking_inequality"], high["endpoint_blocking_inequality"])

    def test_all_Q6_memory_reaches_corner(self):
        memory = corner_target()["physical_memory"]
        self.assertIn("actual_axial_velocity_gd_minus_A1gd", memory)
        self.assertIn("whole_dZ_Jz_plus_c_rho_P", memory)
        self.assertIn("full_pressure_cutoff_derivatives", memory)

    def test_invalid_orders_are_rejected(self):
        for value in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                corner_target(value, 0)
            with self.assertRaises(ValueError):
                endpoint_request(value, 0, 0)

    def test_scope_keeps_endpoint_premise_explicit(self):
        scope = self.result["scope"]
        self.assertTrue(scope["exact_Q1_Q6_residual_is_cartesian_corner_flat_conditionally"])
        self.assertFalse(scope["endpoint_uniform_all_order_physical_bound_proved_here"])
        self.assertFalse(scope["terminal_velocity_or_pressure_regular_at_origin_proved"])
        self.assertFalse(scope["localization_force_energy_comparison_replayed_here"])

    def test_dependency_closure_and_controls(self):
        self.assertEqual(len(self.result["checks"]), 10)
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 6)
        self.assertTrue(all(self.record["gates"].values()))

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["review"]["scope"]["terminal_velocity_or_pressure_regular_at_origin_proved"] = True
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
