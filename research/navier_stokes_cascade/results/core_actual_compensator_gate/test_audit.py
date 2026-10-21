#!/usr/bin/env python3
"""Tests for the actual-compensator data gate and conditional estimates."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import review


class ActualCompensatorGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def test_analytical_profile_contract_is_not_a_formula(self):
        self.assertTrue(all(self.result["parts"]["profile_contract"]["checks"].values()))

    def test_signed_family_preserves_mass_and_second_moment(self):
        rows=self.result["parts"]["signed_same_moment_family"]["checks"]
        self.assertTrue(all(v for k,v in rows.items() if "mass" in k or "second_moment" in k))

    def test_signed_family_shape_primitive_is_unbounded(self):
        rows=self.result["parts"]["signed_same_moment_family"]["checks"]
        self.assertTrue(rows["primitive_has_zero_endpoints"] and rows["primitive_midpoint_nonzero"])
        self.assertTrue(all(v for k,v in rows.items() if "linear_shape_primitive_growth" in k))

    def test_positivity_still_does_not_bound_profile_height(self):
        self.assertTrue(all(self.result["parts"]["positive_concentration"]["checks"].values()))

    def test_exact_slow_composition_constants(self):
        self.assertTrue(all(self.result["parts"]["dilation_budget"]["checks"].values()))
        self.assertEqual(self.result["parts"]["dilation_budget"]["bounds"]["D_Z_B1"],"64")

    def test_parent_pressure_shape_and_radial_results_retained(self):
        self.assertTrue(all(self.result["parts"]["inherited_shape"]["checks"].values()))

    def test_historical_records_and_sources_match(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),3)

    def test_no_scope_promotion(self):
        self.assertFalse(self.record["flags"]["actual_profile_formula_or_numeric_seminorms_found"])
        self.assertFalse(self.record["flags"]["diagnostic_pressure_bump_promoted_to_canonical_data"])
        self.assertFalse(self.record["flags"]["actual_high_target_jet_budget_closed"])
        self.assertFalse(self.record["flags"]["complete_residual_improvement_or_blowup_proved"])

    def test_scope_tampering_is_rejected(self):
        record=deepcopy(self.record);record["flags"]["actual_high_target_jet_budget_closed"]=True
        with self.assertRaises(ValueError):audit.validate(record)

    def test_frozen_record_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"evidence.json";path.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(path)])
            self.assertEqual(path.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
