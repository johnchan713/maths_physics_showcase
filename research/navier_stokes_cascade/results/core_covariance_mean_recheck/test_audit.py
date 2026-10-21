#!/usr/bin/env python3
"""Meaningful exact operator, projection and evidence negative controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import dot, project, review


class MeanDefectRecheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=review()
        cls.record=audit.generate()

    def test_full_radial_operator_against_integrated_pressure(self):
        rows=self.result["parts"]["dynamic"]["checks"]
        for key in ("complete_operator_remainder","plain_pressure_identity","weighted_pressure_identity"):
            self.assertTrue(rows[key],key)

    def test_complete_axial_and_angular_defects(self):
        rows=self.result["parts"]["dynamic"]["checks"]
        self.assertTrue(rows["angular_nonlinear_identity"])
        self.assertTrue(rows["axial_nonlinear_identity"])

    def test_dynamic_omissions_detected(self):
        rows=self.result["parts"]["dynamic"]["checks"]
        for key in ("time_omission_detected","axial_viscosity_omission_detected","axial_transport_omission_detected"):
            self.assertTrue(rows[key],key)

    def test_radial_viscosity_cancels_only_in_both_moments(self):
        rows=self.result["parts"]["dynamic"]["checks"]
        for key in ("radial_vector_viscosity_plain_cancels","radial_vector_viscosity_weighted_cancels","radial_vector_viscosity_pointwise_nonzero"):
            self.assertTrue(rows[key],key)

    def test_later_nonzero_old_mean_is_negative_control(self):
        self.assertTrue(self.result["parts"]["dynamic"]["checks"]["later_nonzero_old_mean_breaks_initial_formula"])

    def test_Haar_products_and_squares_distinguished(self):
        rows=self.result["parts"]["dynamic"]["checks"]
        self.assertTrue(rows["Haar_cross_cancels_only_after_average"])
        self.assertTrue(rows["Haar_square_retained"])

    def test_manufactured_divergence_and_zero_moments(self):
        rows=self.result["parts"]["dynamic"]["checks"]
        for key in ("exact_slow_divergence","both_velocity_moments_zero","compact_endpoint_second_jets"):
            self.assertTrue(rows[key],key)

    def test_projected_bound_saturation_and_strict_unprojected_gap(self):
        self.assertTrue(all(self.result["parts"]["projection"]["checks"].values()))

    def test_exact_small_axial_denominators_grow_energy(self):
        rows=self.result["parts"]["projection"]["small_denominators"]
        for a,b in zip(rows,rows[1:]):
            self.assertGreater(Q(a["D_z"]),Q(b["D_z"]))
            self.assertLess(Q(a["energy"]),Q(b["energy"]))

    def test_zero_tilt_obstruction(self):
        p,d=project([Q(2),Q(2)],[1,1],[1,3])
        self.assertEqual(p,[0,0])
        self.assertEqual(d,0)
        self.assertEqual(dot([3,-1],[2,2],[1,3]),0)

    def test_joint_eliminated_budget_and_dynamic_cancellation_control(self):
        self.assertTrue(all(self.result["parts"]["joint_budget"]["checks"].values()))

    def test_continuum_quarter_tilt_moments(self):
        row=self.result["parts"]["continuum"]
        self.assertEqual(Q(row["D_z"]),Q(99,85))
        self.assertGreater(Q(row["D_P"]),0)

    def test_projection_rejects_inexact_or_invalid_inputs(self):
        with self.assertRaises(TypeError):
            project([1.0,2],[1,1],[1,1])
        for args in (([1],[1],[0]),([1,2],[1],[1,1]),([],[],[]),([1],[0],[1])):
            with self.assertRaises(ValueError):
                project(*args)

    def test_complete_covariance_axial_entries_and_pressure(self):
        rows=self.result["parts"]["complete_covariance"]["checks"]
        for key in ("all_six_covariances_nonzero","omitted_axial_velocity_detected","omitted_axial_pressure_detected","whole_moving_bump_derivative"):
            self.assertTrue(rows[key],key)

    def test_both_physical_pressure_tail_terms_retained(self):
        rows=self.result["parts"]["pressure_tails"]["checks"]
        self.assertTrue(all(rows.values()))
        self.assertTrue(rows["constant_second_profile_can_have_axial_gradient"])
        self.assertTrue(rows["first_tail_cutoff_transition_not_zero"])

    def test_entire_moving_compensator_product_retained(self):
        rows=self.result["parts"]["dynamic"]["checks"]
        for key in ("compensator_unit_mass","full_pressure_primitive_compact","whole_axial_pressure_moment","moving_compensator_derivative_omission_detected"):
            self.assertTrue(rows[key],key)

    def test_all_controls_and_provenance(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)

    def test_no_scientific_scope_promotion(self):
        for key in ("actual_dynamic_terms_cannot_cancel_energy","bounds_uniform_across_all_actual_physical_scales","later_nonzero_means_or_changing_wave_covered_by_initial_formula","all_generated_sources_admissible","complete_residual_improvement_or_blowup_proved","manufactured_tests_are_continuum_proof"):
            self.assertFalse(self.record["flags"][key],key)

    def test_tampered_controls_and_scope_rejected(self):
        for section,key in (("flags","complete_residual_improvement_or_blowup_proved"),("gates","all_exact_mean_controls_pass")):
            altered=deepcopy(self.record)
            altered[section][key]=not altered[section][key]
            with self.assertRaises(ValueError):
                audit.validate(altered)

    def test_frozen_evidence_reproduces(self):
        path=Path(__file__).with_name("evidence.json")
        if path.exists():
            self.assertEqual(json.loads(path.read_text()),self.record)
        else:
            self.fail("Freeze evidence before running the suite")

    def test_output_overwrite_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"record.json"
            path.write_text("sentinel")
            with self.assertRaises(SystemExit):
                audit.main(["--output",str(path)])
            self.assertEqual(path.read_text(),"sentinel")


if __name__=="__main__":
    unittest.main()
