import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q187


class Q187Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q187.statuses().values()))

    def test_actual_symbolic_inputs(self):
        row = q187.actual_symbolic_inputs()
        self.assertTrue(row["all_value_inputs_strictly_positive_and_finite"])
        self.assertEqual(row["u_squared"], "2/epsilon = 2^19*G^22")
        self.assertEqual(row["gamma"], "epsilon/4 = 1/(2^20*G^22)")
        self.assertTrue(row["same_selected_profile"])
        self.assertFalse(row["numerically_evaluated"])

    def test_Cstar_keeps_pressure_and_axial_viscosity(self):
        row = q187.coefficient_contract()
        self.assertTrue(row["actual_shared_budget_and_first_seminorm_functional"])
        self.assertTrue(row["explicit_finite_Cstar"])
        self.assertTrue(row["full_pressure_moment_memory"])
        self.assertTrue(row["pressure_and_axial_viscosity_terms"])
        self.assertFalse(row["numerically_evaluated"])

    def test_exact_selector(self):
        for Cval, rho, h, ell_min in (
            (100, Q(1, 37), Q(1, 200), 2),
            (Q(1, 3), Q(1, 2), Q(1, 17), 99),
            (2**80 + 1, Q(1, 2**30), Q(3, 997), 7),
        ):
            row = q187.exact_dyadic_selector(Cval, rho, h, ell_min)
            self.assertTrue(row["gate"])
            self.assertTrue(row["minimal_dyadic_power"])
            self.assertTrue(row["chart_gate"])
            self.assertGreaterEqual(row["ell"], ell_min)

    def test_selector_contract(self):
        row = q187.selector_contract()
        self.assertTrue(row["rho_strictly_positive_and_finite"])
        self.assertTrue(row["Cstar_strictly_positive_and_finite"])
        self.assertTrue(row["background_certificate_has_exact_radius_formula"])
        self.assertTrue(row["finite_later_selector_exists"])
        self.assertFalse(row["historical_first_integer_or_profile_constant_evaluated"])
        self.assertEqual(
            row["blocking_inequality"],
            "Cstar*2^(-2h*(ell-2)) <= rho_actual",
        )

    def test_frontier_is_primitive_provider_validation(self):
        row = q187.frontier_contract()
        self.assertTrue(row["Q186_actual_symbolic_join_completed"])
        self.assertTrue(row["Q2_source_audited_relative_to_displayed_providers"])
        self.assertTrue(row["modulation_repair_same_map_and_pressure_audited"])
        self.assertFalse(row["primitive_continuation_loop_repair_regional_suprema_independently_validated_here"])
        self.assertFalse(row["numeric_evaluation_required_for_later_tail_existence"])
        self.assertTrue(row["numeric_evaluation_required_to_identify_a_prescribed_historical_band"])

    def test_review_scope(self):
        row = q187.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["actual_symbolic_Q2_selector_closed_relative_to_named_providers"])
        self.assertFalse(row["scope"]["primitive_provider_theorems_reproved_from_initial_equations"])
        self.assertFalse(row["scope"]["uniform_actual_pulse_or_complete_residual_proved"])
        self.assertFalse(row["scope"]["unconditional_force_exclusion_or_blowup_proved"])

    def test_invalid_inputs(self):
        for call in (
            q187.statuses,
            q187.actual_symbolic_inputs,
            q187.coefficient_contract,
            q187.selector_contract,
            q187.frontier_contract,
            q187.review,
            q187.record,
        ):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        for args in ((0, 1, 1), (1, 0, 1), (1, 1, 0)):
            with self.assertRaises(ValueError):
                q187.exact_dyadic_selector(*args)
        with self.assertRaises(TypeError):
            q187.exact_dyadic_selector(1.0, 1, 1)
        with self.assertRaises(TypeError):
            q187.exact_dyadic_selector(1, 1, 1, True)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 433)
        self.assertGreater(record["closure"]["edges"], 1500)
        self.assertGreater(record["closure"]["source_bindings"], 5800)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
