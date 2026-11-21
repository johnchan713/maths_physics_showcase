import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q202


class Q202Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q202.statuses().values()))

    def test_cycle_rates(self):
        r = q202.rates()
        self.assertEqual(r["current_source_D"], Q(74999, 50000))
        self.assertEqual(r["current_target_C0"], Q(99999, 50000))
        self.assertEqual(r["cycle_gain_delta"], Q(49997, 100000))
        self.assertEqual(r["next_source_B1"], Q(39999, 20000))
        self.assertEqual(r["next_required_target_C1"], Q(49999, 20000))
        self.assertEqual(r["proved_target_bridge"], Q(249997, 100000))
        self.assertEqual(r["target_slack"], Q(1, 50000))

    def test_complete_current_source_is_D(self):
        r = q202.current_source()
        self.assertEqual(r["order"], Q(74999, 50000))
        self.assertEqual(r["angular_harmonics"], (-4, -3, -2, -1, 1, 2, 3, 4))
        self.assertTrue(r["full_pressure_gradient_viscosity_cylindrical_crosses_and_squares_retained"])
        self.assertFalse(r["flat_Gaussian_and_base_exterior_terms_forward_inverted"])

    def test_sourced_inverse_orders_and_constraint(self):
        r = q202.sourced_inverse()
        self.assertEqual(r["complete_cut_velocity"], Q(74999, 50000))
        self.assertEqual(r["pressure_and_vector_potential"], Q(99999, 50000))
        self.assertEqual(r["complete_curl_remainder"], Q(199997, 100000))
        self.assertTrue(r["pressure_equals_current_C0"])
        self.assertTrue(r["moving_normal_term_retained"])
        self.assertTrue(r["source_normal_term_retained"])
        self.assertTrue(r["potential_and_pressure_cut_before_complete_curl"])
        self.assertFalse(r["uncut_exit_declared_zero"])

    def test_physical_request_counts_and_domains(self):
        r = q202.physical_domain()
        self.assertEqual(r["A0_integrations"], 30100002)
        self.assertEqual(r["A1_integrations"], 30050003)
        self.assertEqual(r["A0_source_order"], 60200007)
        self.assertEqual(r["A1_source_order"], 60100009)
        self.assertEqual(r["maxima"]["new_covariance"], 60200020)
        self.assertEqual(r["maxima"]["next_signed_target_input"], 60200021)
        self.assertEqual(r["common_wave_pole"], 361200138)
        self.assertEqual(r["post_signed_and_fast_domain"], 722400276)
        self.assertEqual(r["full_nonzero_source_domain"], 1083600414)
        self.assertEqual(r["final_radial_source_domain"], 1444800552)

    def test_physical_compact_margins_are_strict(self):
        r = q202.physical_domain()
        self.assertEqual(r["strict_A0_exponent_margin"], Q(1, 20000000))
        self.assertEqual(r["strict_A1_exponent_margin"], Q(1, 20000000))
        self.assertFalse(r["actual_high_jet_values_evaluated"])

    def test_normalized_mean_domain_is_not_physical_domain(self):
        r = q202.physical_domain()
        self.assertEqual(r["normalized_mean_pole_for_default_request"], 1200114)
        self.assertGreater(r["common_wave_pole"], r["normalized_mean_pole_for_default_request"])
        self.assertFalse(r["normalized_mean_pole_suffices_for_physical_request"])
        self.assertTrue(q202.correction()["physical_domain_strictly_larger"])

    def test_next_coupled_cycle(self):
        r = q202.next_cycle()
        self.assertEqual(set(r["source_limiters"]), {"sourced_full_linear", "signed_full_linear"})
        self.assertEqual(r["next_source"], Q(39999, 20000))
        self.assertEqual(r["next_required_target"], Q(49999, 20000))
        self.assertEqual(r["target_slack"], Q(1, 50000))
        self.assertTrue(r["source_and_target_advanced_simultaneously"])
        self.assertTrue(r["closes_relative_to_named_continuum_inputs"])

    def test_physical_powers_do_not_claim_decay(self):
        r = q202.physical_powers()
        self.assertEqual(r["sourced_velocity_zero_jet"], Q(-4975001, 10000000))
        self.assertEqual(r["sourced_pressure_zero_jet"], Q(-10000001, 10000000))
        self.assertEqual(r["sourced_potential_zero_jet"], Q(49999, 10000000))
        self.assertEqual(r["next_regular_source_zero_jet"], Q(-6000001, 4000000))
        self.assertLess(r["next_regular_source_zero_jet"], 0)

    def test_review_scope(self):
        r = q202.review()
        self.assertTrue(all(r["checks"].values()))
        self.assertTrue(r["scope"]["conditional_analytical_next_sourced_cycle"])
        self.assertFalse(r["scope"]["tests_or_integer_request_arithmetic_are_continuum_proof"])
        self.assertFalse(r["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q202.statuses, q202.current_source, q202.sourced_inverse, q202.next_cycle, q202.correction, q202.frontier, q202.review, q202.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError): q202.rates(0)
        with self.assertRaises(ValueError): q202.physical_domain(-1)
        with self.assertRaises(ValueError): q202.physical_domain(original_pole=0)
        with self.assertRaises(ValueError): q202.physical_powers(Q(1, 50))

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 451)
        self.assertGreaterEqual(record["closure"]["edges"], 1621)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6079)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
