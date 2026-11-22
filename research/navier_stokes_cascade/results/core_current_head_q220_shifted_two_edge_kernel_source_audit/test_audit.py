import json
import tempfile
import unittest
from pathlib import Path
import audit
import q220


class Q220Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q220.statuses().values()))

    def test_directed_maps(self):
        row = q220.directed_maps()
        self.assertTrue(row["physical_shift_retained"])
        self.assertTrue(row["differentiate_at_fixed_integration_shift"])
        self.assertFalse(row["translation_depends_on_differentiated_variables"])
        self.assertEqual(row["extra_M_power_from_differentiation"], 0)
        self.assertFalse(row["wrong_edge_orientation_valid"])

    def test_every_fixed_pole(self):
        for pole in (0, 1, 12, 2_400_066, 2_626_513_221_816):
            row = q220.weight_contract(pole)
            self.assertEqual(row["output_pole"], pole)
            self.assertEqual(row["pole_growth"], 0)
            self.assertTrue(row["all_three_branches_finite_for_fixed_positive_data"])
            self.assertFalse(row["constant_uniform_in_pole"])
        for bad in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                q220.weight_contract(bad)

    def test_cutoff_contract(self):
        row = q220.cutoff_contract()
        self.assertTrue(row["I_minus_selected_away_from_right_edge"])
        self.assertTrue(row["I_plus_selected_away_from_left_edge"])
        self.assertTrue(row["chi_U_support_has_positive_weight_floor"])
        self.assertTrue(row["A_e_keeps_whole_shifted_moment"])
        self.assertFalse(row["zero_Haar_alone_deletes_A_e"])

    def test_typed_boundary(self):
        row = q220.typed_boundary()
        self.assertFalse(row["half_weight_input_automatically_upgraded"])
        self.assertFalse(row["unweighted_jets_plus_support_and_moments_imply_full_weight"])
        self.assertTrue(row["physical_covariance_products_supply_full_weight_by_structure"])
        self.assertTrue(row["source_membership_is_separate_from_operator_boundedness"])

    def test_source_classification(self):
        row = q220.source_classification()
        self.assertTrue(row["Gaussian_absorbs_every_fixed_polynomial_pole"])
        self.assertTrue(row["symbolic_exponential_sum_is_exact_and_finite"])
        self.assertFalse(row["floating_point_materialization_required"])
        self.assertFalse(row["tests_used_as_continuum_proof"])
        self.assertFalse(row["actual_generated_source_membership_proved_by_kernel_itself"])

    def test_review(self):
        row = q220.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q216_two_edge_kernel_quantifier_supported"])
        self.assertFalse(row["scope"]["blanket_full_weight_generated_source_claim"])
        self.assertFalse(row["scope"]["unconditional_force_exclusion_or_blowup"])

    def test_identity_guards(self):
        for call in (q220.statuses, q220.directed_maps, q220.cutoff_contract,
                     q220.typed_boundary, q220.source_classification,
                     q220.review, q220.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q220.weight_contract(1, identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 470)
        self.assertGreaterEqual(record["closure"]["edges"], 1779)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6318)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
