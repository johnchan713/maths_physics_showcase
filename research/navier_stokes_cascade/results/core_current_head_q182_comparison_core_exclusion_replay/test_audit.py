import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q182


class Q182Tests(unittest.TestCase):
    def test_provider_statuses_and_same_object_chain(self):
        self.assertTrue(all(q182.statuses().values()))
        self.assertTrue(all(q182.identity_chain().values()))

    def test_pressure_force_and_harmonic_gradient(self):
        row = q182.pressure_comparison_contract()
        self.assertIn("-div f", row["individual_pressures"])
        self.assertEqual(row["force_cancellations"], 1)
        self.assertTrue(row["finite_energy_essential"])
        self.assertFalse(row["Poisson_alone_removes_harmonic_gradient"])
        self.assertFalse(row["competitor_pressure_decay_assumed"])
        self.assertIn("v=t e1", row["finite_energy_counterexample"])

    def test_commutator_and_radius_ledger(self):
        row = q182.pressure_comparison_contract()
        self.assertEqual(row["commutator_exponent"], Q(4, 3))
        self.assertEqual(row["commutator_radial_coefficient"], 4)
        self.assertEqual(row["commutator_norm_decay"], Q(3, 4))
        self.assertEqual(row["cutoff_power"], 8)
        self.assertFalse(row["seventh_cutoff_power_sufficient"])
        self.assertEqual(row["compact_mixed_flux_decay"], Q(5, 2))
        self.assertEqual(row["young_decays"], (4, Q(14, 5)))
        self.assertEqual(row["worst_radius_decay"], 2)

    def test_comparison_adds_no_competitor_regularities(self):
        row = q182.pressure_comparison_contract()
        self.assertIn("R^-2", row["localized_inequality"])
        self.assertIn("v=u", row["conclusion"])
        self.assertFalse(row["competitor_derivative_bound_assumed"])
        self.assertFalse(row["competitor_dissipation_bound_assumed"])
        self.assertFalse(row["competitor_boundary_flux_bound_assumed"])
        self.assertFalse(row["terminal_uniform_constants_assumed"])

    def test_protected_core_path_and_threshold(self):
        row = q182.protected_core_contract()
        self.assertEqual(row["A"], Q(101, 200))
        self.assertEqual(row["X_in"], Q(1, 101))
        self.assertEqual(row["X_a"], Q(4, 101))
        self.assertEqual(row["support_hole_margin"], Q(3, 101))
        self.assertFalse(row["background_terms_vanish_exactly"])
        self.assertTrue(row["annular_diagonal_terms_vanish_exactly"])
        self.assertTrue(row["threshold_strictly_positive"])
        self.assertTrue(row["normalized_error_within_half"])

    def test_legal_limit_order_and_scope(self):
        row = q182.protected_core_contract()
        self.assertTrue(row["strict_time_order"])
        self.assertIn("R to infinity", row["second_limit"])
        self.assertIn("tau", row["third_limit"])
        self.assertFalse(row["uniform_comparison_constant_as_T_to_one_required"])
        self.assertIn("no globally smooth", row["conditional_conclusion"])
        self.assertFalse(row["future_velocity_or_pressure_constructed"])
        self.assertFalse(row["unconditional_blowup_claimed"])

    def test_review_corrections_and_limits(self):
        row = q182.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["correction"]["Poisson_equation_alone_does_not_remove_harmonic_gradient"])
        self.assertTrue(row["scope"]["P7_same_force_competitor_exclusion_closed_conditionally"])
        self.assertFalse(row["scope"]["all_named_continuum_providers_independently_reproved"])
        self.assertFalse(row["scope"]["unconditional_admissible_force_or_blowup"])

    def test_invalid_identity(self):
        for fn in (
            q182.statuses, q182.identity_chain, q182.pressure_comparison_contract,
            q182.protected_core_contract, q182.review,
        ):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 427)
        self.assertGreater(record["closure"]["edges"], 1476)
        self.assertGreater(record["closure"]["source_bindings"], 5754)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
