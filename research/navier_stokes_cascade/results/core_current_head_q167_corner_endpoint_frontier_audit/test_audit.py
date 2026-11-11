import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q167


class Q167Tests(unittest.TestCase):
    def test_complete_same_object_identity_chain(self):
        row = q167.provider_contracts()
        self.assertTrue(row["all_identity_links_hold"])
        self.assertTrue(all(row["identity_chain"].values()))
        self.assertTrue(row["Q162_same_physical_field_cutoff_support_and_pressure_gauge"])
        self.assertTrue(row["Q154_supplies_Q7_endpoint_blocker_conditionally"])
        self.assertFalse(row["singular_origin_regularized"])

    def test_Q166_residual_reaches_closed_cartesian_corner(self):
        row = q167.corner_contract()
        self.assertEqual(row["selected_h"], Q(1, 200))
        self.assertEqual(row["selected_kappa"], Q(1, 100000))
        self.assertEqual(row["jacobian_floor"], Q(99, 100))
        self.assertEqual(row["sample_stages"], [12817, 18418, 37619])
        self.assertTrue(row["all_controls_hold"])
        self.assertTrue(row["all_six_residual_groups_reach_corner"])
        self.assertTrue(row["complete_memory_reaches_corner"])
        self.assertFalse(row["physical_support_distance_or_edge_weight_floor_used"])

    def test_Q154_supplies_full_endpoint_inequality(self):
        row = q167.endpoint_contract()
        self.assertTrue(row["endpoint_inequality_supplied_conditionally"])
        self.assertEqual(row["implicit_floor"], Q(99, 100))
        self.assertEqual(row["full_closed_band_window"], list(range(18, 27)))
        self.assertTrue(row["observation_inside_qstar"])
        self.assertTrue(row["observation_beyond_first_band"])
        self.assertTrue(row["finite_active_band_box_and_stage_sets_on_each_punctured_compact"])

    def test_endpoint_derivative_ledger_is_exact(self):
        row = q167.endpoint_contract()["localization_request"]
        self.assertEqual(row, {
            "m": 5, "potential": 9, "direct_swirl": 8,
            "pressure": 7, "corner_residual": 6,
        })

    def test_full_exterior_pressure_convolution_is_retained(self):
        row = q167.endpoint_contract()
        self.assertTrue(row["full_pressure_Leibniz_convolution_retained"])
        self.assertEqual(row["pressure_S2_exceeds_two_A2_by"], Q(2) * Q(1, 200)**2 * Q(201, 200)**2)

    def test_zero_Taylor_radius_and_singular_origin_limits(self):
        row = q167.endpoint_contract()
        self.assertTrue(row["Taylor_ratios_grow_unbounded"])
        self.assertFalse(row["exterior_time_Taylor_radius_positive"])
        self.assertTrue(row["compatible_one_sided_fixed_order_limits"])
        self.assertFalse(row["singular_origin_included"])

    def test_frontier_moves_to_Q8_without_localization_promotion(self):
        gate = q167.remaining_gate()
        self.assertIn("closed-conditionally", gate["U4_corner_and_punctured_endpoint"])
        self.assertIn("Q8", gate["next"])
        self.assertFalse(gate["localization_force_energy_comparison_replayed_here"])
        self.assertFalse(gate["unconditional_theorem_gate_closed"])
        self.assertEqual(q167.frontier()[-1], {"id": "U5", "status": "next-Q8-localization-force-audit"})

    def test_scope_rejects_unconditional_promotion(self):
        row = q167.record()
        self.assertTrue(row["U4_corner_endpoint_closed_conditionally"])
        self.assertFalse(row["singular_origin_regularized"])
        self.assertFalse(row["localization_force_energy_comparison_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_identity_rejected(self):
        for fn in (q167.provider_contracts, q167.corner_contract,
                   q167.endpoint_contract, q167.remaining_gate, q167.frontier):
            with self.assertRaises(ValueError):
                fn("wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 406)
        self.assertGreater(record["closure"]["edges"], 1355)
        self.assertGreater(record["closure"]["source_bindings"], 5515)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
