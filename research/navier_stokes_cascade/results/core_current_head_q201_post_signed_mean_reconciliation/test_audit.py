import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q201


class Q201Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q201.statuses().values()))

    def test_corrected_rates(self):
        r = q201.rates()
        self.assertEqual(r["B"], Q(99999, 100000))
        self.assertEqual(r["C"], Q(149999, 100000))
        self.assertEqual(r["D"], Q(74999, 50000))
        self.assertEqual(r["C0"], Q(99999, 50000))

    def test_actual_fast_and_slow_orders(self):
        o = q201.mean_orders()
        self.assertEqual(o["post_signed_Qcirc_rows"], Q(74999, 50000))
        self.assertEqual(o["fast_tangential"], Q(74999, 50000))
        self.assertEqual(o["fast_radial"], Q(124999, 50000))
        self.assertEqual(o["postfast_scalars"], Q(149999, 100000))
        self.assertEqual(o["slow_tangential"], Q(149999, 100000))
        self.assertEqual(o["slow_radial"], Q(249999, 100000))

    def test_regular_mean_rows_close_at_C0(self):
        o = q201.mean_orders()
        self.assertEqual(o["final_tangential"], Q(99999, 50000))
        self.assertEqual(o["final_radial"], Q(99999, 50000))
        self.assertEqual(o["final_regular_mean"], Q(99999, 50000))

    def test_full_difference_margins(self):
        o = q201.mean_orders()
        self.assertEqual(o["tangential_margins_over_C0"]["fast_full_difference"], Q(24999, 50000))
        self.assertEqual(o["tangential_margins_over_C0"]["slow_full_difference"], Q(49999, 100000))

    def test_compact_counts_have_strict_kappa_margin(self):
        o = q201.mean_orders()
        self.assertEqual(o["A1_integrations"], 50001)
        self.assertEqual(o["A0_integrations"], 100000)
        self.assertEqual(o["tangential_margins_over_C0"]["A1_compact_remainder"], Q(1, 100000))
        self.assertEqual(o["radial_margins_over_C0"]["A0_compact_remainder"], Q(1, 100000))

    def test_complete_mean_memory(self):
        m = q201.exact_mean_memory()
        for key in (
            "A1_sign_positive", "fast_all_components_zero_Haar", "fast_divergence_exact",
            "both_velocity_moments_preserved", "slow_uses_current_postfast_scalars",
            "slow_uses_same_fixed_five_map", "whole_reserved_potential_hull_G_zero_required",
            "all_six_covariances_retained", "axial_velocity_retained", "compact_pressure_retained",
            "moving_rho_and_whole_c_rho_P_derivative_retained",
            "viscosity_time_cylindrical_cross_and_square_terms_retained",
        ):
            self.assertTrue(m[key])

    def test_finite_request_graph(self):
        r = q201.compact_requests()
        self.assertEqual(r["A1_source_order"], 100005)
        self.assertEqual(r["A0_source_order"], 200003)
        self.assertEqual(r["maxima"]["post_signed_row"], 200013)
        self.assertEqual(r["maxima"]["new_covariance"], 200016)
        self.assertEqual(r["maxima"]["next_signed_target_input"], 200017)
        self.assertEqual(r["common_wave_pole"], 1200114)
        self.assertFalse(r["actual_high_jets_evaluated"])

    def test_precise_conditional_inputs(self):
        r = q201.blocking_inputs()
        self.assertIn("Qcirc gr_final", r["A0_inequality"])
        self.assertIn("Qcirc EpostSigned_z", r["A1_inequality"])
        self.assertTrue(r["closed_relative_to_named_continuum_providers"])
        self.assertFalse(r["independently_evaluated_or_formally_certified"])

    def test_physical_residual_still_not_decay(self):
        self.assertEqual(q201.physical_residual_power(), Q(-15025001, 10000000))
        self.assertLess(q201.physical_residual_power(), 0)

    def test_review_scope(self):
        r = q201.review()
        self.assertTrue(all(r["checks"].values()))
        self.assertTrue(r["scope"]["conditional_analytical_first_post_signed_mean_block"])
        self.assertFalse(r["scope"]["tests_or_finite_request_arithmetic_are_continuum_proof"])
        self.assertFalse(r["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q201.statuses, q201.exact_mean_memory, q201.frontier, q201.review, q201.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError): q201.rates(0)
        with self.assertRaises(ValueError): q201.compact_requests(-1)
        with self.assertRaises(ValueError): q201.compact_requests(original_pole=0)
        with self.assertRaises(ValueError): q201.physical_residual_power(Q(1, 50))

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 450)
        self.assertGreaterEqual(record["closure"]["edges"], 1612)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6066)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
