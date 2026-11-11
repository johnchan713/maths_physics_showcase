import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q169


class Q169Tests(unittest.TestCase):
    def test_same_object_identity_chain(self):
        row = q169.provider_contracts()
        self.assertTrue(row["all_identity_links_hold"])
        self.assertTrue(all(row["identity_chain"].values()))
        self.assertTrue(row["same_Q162_field_force_pressure_gauge_and_localization"])
        self.assertTrue(row["Q168_localization_force_energy_closed_conditionally"])

    def test_force_divergence_and_pressure_gauge_are_retained(self):
        row = q169.pressure_contract()
        self.assertIn("-div f", row["individual_pressures"])
        self.assertEqual(row["force_cancellations"], 1)
        self.assertIn("p=c p0", row["constructed_pressure_gauge"])
        self.assertIn("R_i R_j", row["recovery"])
        self.assertFalse(row["competitor_pressure_decay_assumed"])

    def test_harmonic_gradient_is_removed_by_momentum_not_Poisson_alone(self):
        row = q169.pressure_contract()
        self.assertIn("H^-2", row["riesz_pressure_space"])
        self.assertIn("H^-3", row["averaged_gradient_space"])
        self.assertTrue(row["finite_energy_essential"])
        self.assertFalse(row["Poisson_equation_alone_removes_harmonic_gradient"])
        self.assertIn("v=t e1", row["finite_energy_counterexample"])

    def test_exact_commutator_radial_integral(self):
        row = q169.commutator_contract()
        self.assertEqual(row["exponent"], Q(4, 3))
        self.assertEqual(row["near_radial_power"], Q(-2, 3))
        self.assertEqual(row["far_radial_power"], -2)
        self.assertEqual(row["near_coefficient"], 3)
        self.assertEqual(row["far_coefficient"], 1)
        self.assertEqual(row["radial_coefficient"], 4)
        self.assertEqual(row["norm_decay"], Q(3, 4))

    def test_eighth_cutoff_power_is_minimal_for_this_pairing(self):
        row = q169.commutator_contract()
        self.assertTrue(row["cutoff_sufficient"])
        self.assertFalse(row["seventh_power_sufficient"])
        self.assertEqual(row["cutoff_power"], 8)

    def test_compact_mixed_pressure_flux_decay(self):
        row = q169.radius_contract()
        self.assertEqual(row["mixed_pressure_point_decay"], 3)
        self.assertEqual(row["gradient_cutoff_decay"], 1)
        self.assertEqual(row["annular_L1_volume_cost"], Q(3, 2))
        self.assertEqual(row["compact_mixed_flux_decay"], Q(5, 2))

    def test_young_rows_and_limiting_remainder(self):
        row = q169.radius_contract()
        self.assertEqual(
            [item["young_decay"] for item in row["terms"]],
            [4, Q(14, 5)],
        )
        self.assertEqual(
            [item["sobolev_offset_decay"] for item in row["terms"]],
            [Q(5, 2), Q(5, 2)],
        )
        self.assertEqual(row["cutoff_laplacian_decay"], 2)
        self.assertEqual(row["worst_decay"], 2)
        self.assertEqual(row["obsolete_unsplit_decay"], Q(4, 3))

    def test_fixed_T_comparison_adds_no_competitor_regularities(self):
        row = q169.radius_contract()
        self.assertIn("R^-2", row["localized_inequality"])
        self.assertIn("v=u", row["conclusion"])
        self.assertFalse(row["competitor_derivative_bound_assumed"])
        self.assertFalse(row["competitor_dissipation_bound_assumed"])
        self.assertFalse(row["competitor_boundary_flux_bound_assumed"])
        self.assertFalse(row["terminal_uniform_constants_assumed"])

    def test_frontier_moves_only_to_protected_core_audit(self):
        row = q169.remaining_gate()
        self.assertIn("closed-conditionally", row["U6_fixed_gauge_pressure_and_comparison"])
        self.assertIn("protected physical core path", row["next"])
        self.assertFalse(row["protected_core_or_terminal_exclusion_replayed_here"])
        self.assertFalse(row["unconditional_theorem_gate_closed"])
        self.assertEqual(q169.frontier()[-1], {
            "id": "U7", "status": "next-Q157-Q163-protected-core-limit-audit",
        })

    def test_scope_rejects_terminal_promotion(self):
        row = q169.record()
        self.assertTrue(row["U6_pressure_comparison_closed_conditionally"])
        self.assertFalse(row["protected_core_exclusion_or_blowup_certified"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_identity_rejected(self):
        for fn in (q169.provider_contracts, q169.pressure_contract,
                   q169.commutator_contract, q169.radius_contract,
                   q169.remaining_gate, q169.frontier):
            with self.assertRaises(ValueError):
                fn("wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 408)
        self.assertGreater(record["closure"]["edges"], 1364)
        self.assertGreater(record["closure"]["source_bindings"], 5530)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
