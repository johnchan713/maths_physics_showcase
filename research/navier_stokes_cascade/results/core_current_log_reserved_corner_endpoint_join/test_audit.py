from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from transfer import (
    cartesian_constant, closed_chart, corner_target, endpoint_request,
    endpoint, residual, review,
)


class CurrentLogReservedCornerEndpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_closed_chart_exact_geometry(self):
        chart = closed_chart()
        self.assertEqual(chart["terminal_identity"], 1)
        self.assertEqual(chart["jacobian_floor"], Q(99, 100))
        self.assertEqual(chart["normalized_support_margin"], Q(3, 8012))
        self.assertFalse(chart["physical_support_floor_used"])
        self.assertFalse(chart["edge_weight_floor_used"])

    def test_current_stage_and_extra_derivatives_reach_corner(self):
        row = corner_target(1, 3, (7, 11), 13, 17)
        request = row["current_stage_request"]
        self.assertEqual(row["comparison_stage"], request["stage"])
        for name in ("A0", "A1"):
            self.assertEqual(
                request["raw_sources_with_extra_integration"][name],
                request["raw_sources_before_extra_integration"][name] + 2,
            )

    def test_default_corner_uses_current_log_reserved_stage(self):
        row = corner_target()
        parent = residual.comparison()
        self.assertEqual(row["comparison_stage"], 12817)
        self.assertEqual(row["comparison_stage"], parent["log_reserved_comparison_stage"])
        self.assertGreaterEqual(row["closed_chart_exponent"], 1)

    def test_arbitrary_cartesian_powers_close(self):
        for h in (Q(1, 200), Q(9, 1000)):
            for m, power in ((0, 0), (3, 5)):
                row = corner_target(m, power, h=h)
                self.assertTrue(all(row["checks"].values()))
                self.assertGreaterEqual(row["retained_exponent_margin"], 0)
                self.assertEqual(
                    row["cartesian_factor_at_requested_exponent"],
                    4 ** (power + 1),
                )

    def test_all_six_groups_reserves_and_pressure_memory_retained(self):
        row = corner_target(1, 3, (7, 11), 13, 17)
        self.assertEqual(tuple(row["residual_groups"]), residual.RESIDUAL_GROUPS)
        self.assertEqual(set(row["reserve_accounts"]), {"A0", "A1", "regular", "cutoff"})
        self.assertEqual(tuple(row["pressure_memory"]), residual.MEMORY)

    def test_endpoint_default_derivative_ledger(self):
        row = endpoint_request()
        self.assertEqual(row["punctured_field_Cauchy_orders"], {
            "potential": 4, "direct_swirl": 3, "pressure": 2,
        })
        self.assertEqual(
            row["punctured_field_Cauchy_orders"],
            row["localization_Cauchy_orders"],
        )
        self.assertEqual(row["residual_Cauchy_jet"], 1)

    def test_endpoint_current_ledger_matches_every_fixed_request(self):
        for a, j, power in ((0, 0, 0), (1, 3, 5), (5, 4, 11), (3, 7, 2)):
            row = endpoint_request(a, j, power)
            total = a + j
            current = endpoint.requirements(total)
            self.assertEqual(row["residual_Cauchy_jet"], total + 1)
            self.assertEqual(
                current["future_force_Cauchy_potential_bound"], total + 4,
            )
            self.assertTrue(all(row["checks"].values()))
        # The endpoint offsets themselves are exact at arbitrarily high fixed
        # order and do not require evaluating the expensive residual recurrence.
        self.assertEqual(
            endpoint.requirements(200)["future_force_Cauchy_pressure_bound"],
            202,
        )

    def test_request_dependent_stage_does_not_change_field(self):
        low = endpoint_request(0, 0, 0)
        high = endpoint_request(5, 5, 10)
        self.assertGreater(high["comparison_stage"], low["comparison_stage"])
        self.assertEqual(
            low["endpoint_blocking_inequality"],
            high["endpoint_blocking_inequality"],
        )

    def test_invalid_inputs_are_rejected(self):
        for value in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                endpoint_request(value, 0, 0)
        with self.assertRaises(ValueError):
            closed_chart(Q(1, 100))
        with self.assertRaises(ValueError):
            cartesian_constant(-1)

    def test_review_scope_and_controls(self):
        self.assertEqual(len(self.result["checks"]), 10)
        self.assertTrue(all(self.result["checks"].values()))
        scope = self.result["scope"]
        self.assertTrue(scope["current_same_field_residual_cartesian_corner_flat_conditionally"])
        self.assertTrue(scope["current_closed_time_localization_inputs_joined_conditionally"])
        self.assertFalse(scope["endpoint_uniform_all_order_physical_bound_proved_here"])
        self.assertFalse(scope["localization_force_energy_comparison_replayed_here"])

    def test_dependency_provenance_and_tampering(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 6)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))
        self.assertTrue(all(self.record["gates"].values()))
        changed = deepcopy(self.record)
        changed["flags"]["terminal_velocity_or_pressure_regular_at_origin_proved"] = True
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
