import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q219


class Q219Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q219.statuses().values()))

    def test_regional_ledger(self):
        row = q219.regional_ledger()
        self.assertLess(max(row["activation_a_floor_power"],
                            row["activation_third_gap_power"],
                            row["activation_fourth_gap_power"],
                            row["joined_gap_power"]), row["compact_envelope_power"])
        self.assertEqual(row["corrected_outer_ramp_floor"], Q(17701, 30625))
        self.assertTrue(row["corrected_outer_ramp_exceeds_half"])
        self.assertTrue(row["axial_control_exceeds_one"])
        self.assertTrue(row["intermediate_control_exceeds_one"])
        self.assertFalse(row["sample_grid_used_as_continuum_proof"])

    def test_edge_contract(self):
        row = q219.edge_contract()
        self.assertFalse(row["global_absolute_stress_infimum_positive"])
        self.assertEqual(row["global_absolute_stress_infimum"], Q(0))
        self.assertTrue(row["inner_factor_floor_positive"])
        self.assertTrue(row["outer_b0_positive"])
        self.assertTrue(row["normalized_direction_extends_smoothly"])

    def test_background_contract(self):
        row = q219.background_contract()
        self.assertTrue(row["all_fixed_derivative_orders"])
        self.assertTrue(row["one_domain_for_all_fixed_orders"])
        self.assertFalse(row["constants_uniform_in_order"])
        self.assertTrue(row["finite_head_plus_geometric_tail"])
        self.assertTrue(row["coefficient_one_kept_in_finite_head"])
        self.assertFalse(row["generic_absolute_weighted_bound_is_sufficient"])

    def test_perturbation_contract(self):
        row = q219.perturbation_contract()
        self.assertTrue(row["rho_gap_positive_from_finite_positive_inputs"])
        self.assertTrue(row["selector_uses_finitely_many_value_direction_rows"])
        self.assertTrue(row["selector_is_derivative_order_independent"])
        self.assertFalse(row["scale_selection_circularity"])

    def test_source_classification(self):
        row = q219.source_classification()
        self.assertTrue(row["regional_cover_is_continuum_and_finite"])
        self.assertTrue(row["weighted_relative_stress_is_continuum"])
        self.assertFalse(row["diagnostic_only_continuum_row_found"])
        self.assertFalse(row["every_primitive_regional_ODE_rederived_here"])

    def test_review(self):
        row = q219.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q216_same_band_gap_leaf_supported_by_continuum_source_chain"])
        self.assertTrue(row["scope"]["false_positive_absolute_stress_floor_rejected"])
        self.assertFalse(row["scope"]["uniform_actual_pulse_or_complete_residual"])
        self.assertFalse(row["scope"]["unconditional_force_exclusion_or_blowup"])

    def test_identity_guards(self):
        for call in (q219.statuses, q219.regional_ledger, q219.edge_contract,
                     q219.background_contract, q219.perturbation_contract,
                     q219.source_classification, q219.review, q219.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 469)
        self.assertGreaterEqual(record["closure"]["edges"], 1770)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6300)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
