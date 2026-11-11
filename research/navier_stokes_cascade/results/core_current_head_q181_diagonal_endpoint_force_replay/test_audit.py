import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q181


class Q181Tests(unittest.TestCase):
    def test_provider_statuses_and_identity_chain(self):
        self.assertTrue(all(q181.statuses().values()))
        self.assertTrue(all(q181.identity_chain().values()))

    def test_complete_diagonal_residual(self):
        row = q181.diagonal_residual_contract()
        self.assertEqual(row["complete_residual_groups"], (
            "time", "viscosity", "pressure", "old_new", "new_old", "new_square"
        ))
        self.assertTrue(row["all_six_groups_retain_complete_memory"])
        self.assertTrue(row["fixed_jet_target_then_finite_stage"])
        self.assertFalse(row["stage_uniform_constant_or_J_of_q_used"])
        self.assertEqual(row["regular_and_cross_margin"], Q(49997, 40000000))
        self.assertEqual(row["A0_own_log_margin"], Q(1, 40000000))
        self.assertTrue(row["all_fixed_cartesian_residual_jets_flat"])

    def test_endpoint_is_separate_and_nonanalytic(self):
        row = q181.endpoint_contract()
        self.assertEqual(row["jacobian_floor"], Q(99, 100))
        self.assertTrue(row["endpoint_inequality_supplied_conditionally"])
        self.assertTrue(row["full_pressure_Leibniz_convolution_retained"])
        self.assertGreater(row["pressure_S2_excess"], 0)
        self.assertFalse(row["exterior_time_Taylor_radius_positive"])
        self.assertFalse(row["singular_origin_included"])

    def test_exact_localization_force_and_energy(self):
        row = q181.localization_force_energy_contract()
        self.assertEqual(row["term_count"], 12)
        self.assertIn("p0 grad c", row["terms"])
        self.assertIn("(d.grad)d", row["terms"])
        self.assertIn("=0", row["vanishing_thirteenth"])
        self.assertTrue(row["divergence_free"] and row["zero_initial_velocity"])
        self.assertTrue(row["one_compact_smooth_force"])
        self.assertTrue(row["factor_two_retained"])
        self.assertFalse(row["combined_endpoint_extrema_budget_claimed"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])

    def test_review_corrections_and_scope(self):
        row = q181.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["correction"]["residual_flatness_does_not_supply_endpoint_field_bounds"])
        self.assertTrue(row["correction"]["zero_energy_division_is_replaced_by_epsilon_regularization"])
        self.assertTrue(row["scope"]["P6_closed_relative_to_named_continuum_providers"])
        self.assertFalse(row["scope"]["terminal_origin_regularized"])
        self.assertFalse(row["scope"]["global_residual_improvement_proved"])
        self.assertFalse(row["scope"]["comparison_core_exclusion_or_blowup_replayed"])

    def test_invalid_identity(self):
        for fn in (
            q181.statuses, q181.identity_chain, q181.diagonal_residual_contract,
            q181.endpoint_contract, q181.localization_force_energy_contract, q181.review,
        ):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 426)
        self.assertGreater(record["closure"]["edges"], 1466)
        self.assertGreater(record["closure"]["source_bindings"], 5740)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
