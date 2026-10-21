#!/usr/bin/env python3
"""Exact full-shape, positivity, pressure and no-scope-promotion checks."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import expr,eadd,escale,gshift,fixture,review


class PressureShapeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.f=fixture();cls.result=review();cls.record=audit.generate()

    def test_reference_density_mass_and_moment(self):
        for key in ("reference_mass_one","reference_second_moment_gate","reference_gap_density_exact"):
            self.assertTrue(self.f["checks"][key])

    def test_same_moment_perturbation_has_certified_positivity(self):
        for key in ("perturbation_support_inside_gap","perturbation_zero_mass","perturbation_zero_second_moment","exact_positivity_budget"):
            self.assertTrue(self.f["checks"][key])

    def test_pointwise_pressure_survives_vanishing_integrated_coefficient(self):
        for key in ("integrated_fourth_coefficient_zero","pointwise_fourth_coefficient_nonzero","compact_pressure_primitive","pressure_primitive_derivative"):
            self.assertTrue(self.f["checks"][key])

    def test_compact_signed_target_does_not_remove_shape_mismatch(self):
        for key in ("compact_signed_target_primitive","signed_target_derivative","signed_target_pointwise_nonzero"):
            self.assertTrue(self.f["checks"][key])

    def test_complete_radial_operator_before_shape_compensation(self):
        rows={k:v for k,v in self.result["parts"]["pointwise_operator"]["checks"].items() if "complete_source_reference_cancellation" in k}
        self.assertEqual(len(rows),5);self.assertTrue(all(rows.values()))

    def test_actual_compensated_pressure_on_every_piece(self):
        rows={k:v for k,v in self.result["parts"]["pointwise_operator"]["checks"].items() if "pointwise_compensated_pressure" in k}
        self.assertEqual(len(rows),5);self.assertTrue(all(rows.values()))

    def test_nonpressure_axial_velocity_terms_have_lower_target_order(self):
        rows={k:v for k,v in self.result["parts"]["pointwise_operator"]["checks"].items() if "all_nonpressure" in k}
        self.assertEqual(len(rows),5);self.assertTrue(all(rows.values()))

    def test_radial_high_jet_survives_any_unit_compensation(self):
        rows=self.result["parts"]["pointwise_operator"]["checks"]
        self.assertTrue(rows["radial_high_jet_total_nonzero"])
        self.assertTrue(all(v for k,v in rows.items() if "radial_compensation_identity" in k))

    def test_full_pressure_is_compact_but_not_pointwise_zero(self):
        rows=self.result["parts"]["pointwise_operator"]["checks"]
        self.assertTrue(rows["complete_pressure_right_endpoint_zero"])
        self.assertTrue(rows["integrated_axial_high_jet_zero_but_pointwise_nonzero"])

    def test_all_moving_density_product_terms_kept(self):
        rows=self.result["parts"]["moving_density"]["checks"]
        self.assertTrue(rows["full_shape_commutator_identity"])
        self.assertTrue(rows["every_moving_shape_term_nonzero"])
        self.assertTrue(rows["shape_match_leaves_commutator"])

    def test_moving_pressure_radial_and_axial_identities(self):
        rows=self.result["parts"]["moving_density"]["checks"]
        for key in ("full_compensated_source_zero_total","full_pressure_compact","full_pressure_radial_identity","entire_moving_pressure_moment","moving_compensator_product_nonzero"):
            self.assertTrue(rows[key])

    def test_half_L1_primitive_bound_has_a_sharp_gap_control(self):
        self.assertTrue(all(self.result["parts"]["primitive_bound"]["checks"].values()))

    def test_log_normalization_never_needs_a_float(self):
        value=expr(one=2,g=3,b=-7)
        self.assertEqual(gshift(gshift(value,-1),1),value)
        self.assertEqual(eadd(value,escale(value,-1)),{})
        self.assertEqual(self.f["controls"]["reference_mass"],{"logg^0*logb^0":"1"})

    def test_complete_wave_covariance_and_axial_velocity_retained(self):
        rows=self.result["parts"]["complete_covariance"]["checks"]
        self.assertTrue(rows["all_six_covariances_nonzero"])
        self.assertTrue(rows["omitted_axial_velocity_detected"])

    def test_both_full_background_pressure_tails_retained(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))

    def test_historical_records_and_source_hashes_match(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)

    def test_reference_density_does_not_replace_actual_compensator(self):
        for key in ("actual_compensator_or_any_frozen_profile_changed","actual_compensator_shape_match_verified","actual_linked_canonical_sequence_refuted","actual_high_target_jet_bounds_or_full_residual_closed","unconditional_blowup_proved","finite_diagnostics_are_actual_continuum_suprema"):
            self.assertFalse(self.record["flags"][key])

    def test_scientific_scope_tampering_rejected(self):
        record=deepcopy(self.record);record["flags"]["actual_compensator_shape_match_verified"]=True
        with self.assertRaises(ValueError):audit.validate(record)

    def test_frozen_record_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"evidence.json";path.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(path)])
            self.assertEqual(path.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
