import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q199


class Q199Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q199.statuses().values()))

    def test_corrected_current_rates(self):
        row = q199.rates()
        self.assertEqual(row["B"], Q(99999, 100000))
        self.assertEqual(row["C"], Q(149999, 100000))
        self.assertEqual(row["D"], Q(74999, 50000))
        self.assertEqual(row["C0"], Q(99999, 50000))
        self.assertEqual(row["conditional_target_gain"], Q(49999, 100000))

    def test_polynomial_fixture_mass_positive(self):
        self.assertGreater(q199.polynomial_bump_mass(), 0)

    def test_scaled_bumps_are_disjoint_interior(self):
        for eta in (Q(1, 64), Q(1, 32), Q(1, 16)):
            row = q199.bump_geometry(eta)
            self.assertTrue(row["supports_disjoint"])
            self.assertTrue(row["supports_strictly_interior"])
            self.assertEqual(row["signed_total_mass"], 0)
            self.assertGreater(row["primitive_at_observation_without_common_amplitude"], 0)

    def test_every_fixed_pole_has_negative_power(self):
        for pole in range(101):
            row = q199.pole_obstruction(pole)
            self.assertEqual(row["weight_exponential_gap"], Q(7, 36))
            self.assertLess(row["power_after_lower_bound"], 0)
            self.assertTrue(row["defeats_fixed_pole"])

    def test_fixed_unweighted_jets_can_tend_to_zero(self):
        for order in range(101):
            row = q199.pole_obstruction(order)
            self.assertGreater(row["fixed_unweighted_derivative_decay_power"], 0)
            self.assertTrue(row["fixed_unweighted_derivatives_can_tend_to_zero"])

    def test_dependency_cut_has_no_actual_certificate(self):
        row = q199.provider_cut()
        self.assertFalse(row["all_required_actual_certificates_present"])
        self.assertFalse(any(row.values()))

    def test_requested_inequality_and_scope(self):
        request = q199.requested_inequality()
        self.assertEqual(request["components"], (1, 2))
        self.assertEqual(request["rate"], Q(149999, 100000))
        self.assertFalse(request["proved_here"])
        review = q199.review()
        self.assertTrue(all(review["checks"].values()))
        self.assertFalse(review["scope"]["actual_weighted_current_target_or_signed_pulse_certified"])

    def test_physical_limit_unchanged(self):
        self.assertEqual(q199.physical_residual_power(), Q(-15025001, 10000000))
        self.assertLess(q199.physical_residual_power(), 0)

    def test_invalid_inputs(self):
        for call in (q199.statuses, q199.frontier, q199.review, q199.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError): q199.rates(0)
        with self.assertRaises(ValueError): q199.bump_geometry(0)
        with self.assertRaises(ValueError): q199.bump_geometry(Q(1, 2))
        with self.assertRaises(ValueError): q199.pole_obstruction(-1)
        with self.assertRaises(ValueError): q199.physical_residual_power(Q(1, 50))

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 448)
        self.assertGreaterEqual(record["closure"]["edges"], 1602)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6048)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
