#!/usr/bin/env python3
"""Exact potential-gap, complete mean-operator and pressure-jet controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import integrate_laurent,integrate_xlog,logs,ladd,lscale,means,radial_fixture,review


class SlowJetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.radial=radial_fixture();cls.result=review();cls.record=audit.generate()

    def test_prescribed_axial_flux_and_current_target_rows(self):
        for key in ("axial_zero_flux","fractional_moment_dilation","exact_current_angular_target"):
            self.assertTrue(self.radial["checks"][key])

    def test_potential_moment_from_piecewise_primitive_and_source_log_row(self):
        for key in ("logarithmic_potential_moment","independent_log_source_row","weighted_potential_moment"):
            self.assertTrue(self.radial["checks"][key])

    def test_zero_flux_does_not_mean_zero_potential_integral(self):
        self.assertGreater(self.radial["Hlog"],0)
        self.assertTrue(self.radial["checks"]["no_spurious_log_endpoint_or_constant"])

    def test_gap_is_part_of_actual_potential_support(self):
        self.assertTrue(self.radial["checks"]["gap_potential_is_nonzero"])
        self.assertTrue(self.radial["checks"]["omitting_gap_changes_moment"])
        self.assertTrue(self.radial["checks"]["both_potential_endpoints_zero"])

    def test_elementary_laurent_and_log_integrals(self):
        R=means.R
        self.assertEqual(integrate_laurent(R.shift(-2),Q(1),Q(2),logs(),logs(g=1)),logs(g=1))
        self.assertEqual(integrate_xlog(means.Poly(1),Q(1),Q(2),logs(),logs(g=1)),logs(one=-Q(3,4),g=2))
        with self.assertRaises(ValueError):integrate_laurent(R*means.Z,Q(1),Q(2),logs(),logs(g=1))

    def test_odd_part_of_complete_radial_operator_isolates_cubic_response(self):
        checks=self.result["parts"]["full_operator"]["checks"]
        rows={k:v for k,v in checks.items() if "full_operator_principal" in k}
        self.assertEqual(len(rows),6);self.assertTrue(all(rows.values()))

    def test_axial_viscosity_cannot_be_deleted_after_radial_moment_cancellation(self):
        checks=self.result["parts"]["full_operator"]["checks"]
        rows={k:v for k,v in checks.items() if "axial_viscosity_omission" in k}
        self.assertEqual(len(rows),6);self.assertTrue(all(rows.values()))

    def test_full_axial_moment_retains_complete_pressure_response(self):
        checks=self.result["parts"]["full_operator"]["checks"]
        for degree in (3,4):
            self.assertTrue(checks[f"complete_axial_moment_odd_part_degree{degree}"])
            self.assertTrue(checks[f"pressure_principal_coefficient_degree{degree}"])

    def test_Q_normalization_exponents_include_s_to_h(self):
        checks=self.result["parts"]["full_operator"]["checks"]
        rows={k:v for k,v in checks.items() if "moment_exponent" in k}
        self.assertEqual(len(rows),4);self.assertTrue(all(rows.values()))

    def test_special_fourth_coefficient_gate_is_exact_but_not_assumed_actual(self):
        checks=self.result["parts"]["pressure_gate"]["checks"]
        self.assertTrue(checks["exact_fourth_order_cancellation_gate"])
        self.assertTrue(checks["algebraically_tuned_principal_coefficient_zero"])
        self.assertFalse(self.record["flags"]["actual_compensator_special_cancellation_identity_verified"])

    def test_cancelling_fourth_coefficient_does_not_cancel_pressure_third_coefficient(self):
        self.assertTrue(self.result["parts"]["pressure_gate"]["checks"]["pressure_third_coefficient_not_cancelled"])

    def test_whole_moving_pressure_commutator_survives(self):
        checks=self.result["parts"]["pressure_gate"]["checks"]
        self.assertTrue(checks["full_moving_pressure_commutator"])
        self.assertTrue(checks["moving_pressure_commutator_survives"])

    def test_same_order_norm_frequency_control_and_inverse_tilt_denominator(self):
        self.assertTrue(all(self.result["parts"]["rate_control"]["checks"].values()))

    def test_complete_six_entry_covariance_and_axial_velocity_retained(self):
        checks=self.result["parts"]["complete_covariance"]["checks"]
        self.assertTrue(checks["all_six_covariances_nonzero"])
        self.assertTrue(checks["omitted_axial_velocity_detected"])

    def test_both_background_pressure_tails_retained(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))

    def test_symbolic_logs_cancel_without_numerical_roundoff(self):
        value=logs(one=Q(2),g=Q(7),b=Q(-3))
        self.assertEqual(ladd(value,lscale(value,-1)),{})
        self.assertEqual(set(self.radial["controls"]["potential_moment"]),{"g"})

    def test_all_historical_record_and_source_bindings_pass(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)

    def test_isolated_map_obstruction_not_promoted_to_canonical_failure(self):
        for key in ("all_smooth_finite_cycles_or_actual_canonical_sequence_refuted","actual_uniform_high_target_jet_budgets_proved","new_tail_profile_or_tilt_chosen","unconditional_residual_improvement_or_blowup_proved","finite_controls_are_actual_continuum_suprema"):
            self.assertFalse(self.record["flags"][key])

    def test_scientific_scope_tampering_rejected(self):
        record=deepcopy(self.record);record["flags"]["actual_uniform_high_target_jet_budgets_proved"]=True
        with self.assertRaises(ValueError):audit.validate(record)

    def test_frozen_record_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_output_overwrite_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"evidence.json";path.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(path)])
            self.assertEqual(path.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
