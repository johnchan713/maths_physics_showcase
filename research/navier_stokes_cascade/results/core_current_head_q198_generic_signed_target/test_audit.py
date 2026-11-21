import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q198


class Q198Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q198.statuses().values()))

    def test_corrected_rates_and_pure_ledger(self):
        row = q198.rates()
        self.assertEqual(row["B_current_error_and_signed_principal_log"], Q(99999, 100000))
        self.assertEqual(row["C_target_and_covariance_log"], Q(149999, 100000))
        self.assertEqual(row["D_complete_curl_and_current_source_log"], Q(74999, 50000))
        self.assertEqual(row["C0_covariance_error_and_projected_target_log"], Q(99999, 50000))
        self.assertEqual(row["target_gain_log"], Q(49999, 100000))
        self.assertEqual(row["C0_pure"], Q(199997, 100000))

    def test_complete_covariance_orders(self):
        row = q198.covariance_orders()
        self.assertTrue(row["all_errors_tie_at_C0"])
        self.assertEqual(set(row["error_terms"].values()), {Q(99999, 50000)})
        self.assertEqual(row["gain_condition_minimum"], Q(49999, 100000))

    def test_exact_generic_complex_covariance_identity(self):
        row = q198.review()["covariance_fixture"]
        self.assertTrue(row["identity_exact"])
        self.assertTrue(row["all_six_symmetric_entries_present"])
        self.assertTrue(row["complex_quadratures_retained"])
        self.assertTrue(row["unequal_harmonics_not_crossed"])

    def test_pressure_projection_orders_not_promoted(self):
        row = q198.pressure_projection_orders()
        self.assertEqual(row["Haar_pressure_change"], Q(149999, 100000))
        self.assertEqual(row["pointwise_pressure_change"], Q(74999, 50000))
        self.assertEqual(row["projected_tangential_target"], Q(99999, 50000))
        self.assertFalse(row["pressure_or_scalar_promoted_to_C0"])

    def test_moving_projection_and_pressure_terms(self):
        row = q198.moving_derivatives(2, 3, 5, 7, 11, 13, 17, 19, 23)
        self.assertEqual(row["projection_commutator"], 30)
        self.assertEqual(row["whole_pressure_product_derivative"], 658)
        self.assertTrue(row["moving_terms_nonzero"])

    def test_enlarged_pressure_domain(self):
        row = q198.pressure_domain(0, 1, 12)
        self.assertEqual(row["common_wave_pole"], 30)
        self.assertEqual(row["covariance_pressure_source_pole"], 60)
        self.assertFalse(row["old_domain_suffices"])

    def test_unweighted_bound_does_not_supply_edge_weight(self):
        small = q198.edge_counterexample(10)
        large = q198.edge_counterexample(100)
        self.assertEqual(small["unweighted_value"], 1)
        self.assertGreater(large["weighted_quotient"], small["weighted_quotient"])

    def test_scope_and_physical_limit(self):
        row = q198.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertEqual(row["physical_residual_power"], Q(-15025001, 10000000))
        self.assertFalse(row["frontier"]["actual_weighted_input_inequality_proved_here"])
        self.assertFalse(row["scope"]["actual_signed_pulse_or_post_signed_mean_refresh_certified"])

    def test_invalid_inputs(self):
        for call in (q198.statuses, q198.frontier, q198.review, q198.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError): q198.rates(0)
        with self.assertRaises(ValueError): q198.pressure_domain(0, 0, 1)
        with self.assertRaises(ValueError): q198.pressure_domain(-1, 1, 1)
        with self.assertRaises(ValueError): q198.edge_counterexample(1)
        with self.assertRaises(ValueError): q198.physical_residual_power(Q(1, 50))
        with self.assertRaises(ValueError): q198.covariance({0: ((1, 2, 3), (1, 2, 3))})

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 447)
        self.assertGreaterEqual(record["closure"]["edges"], 1598)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6040)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
