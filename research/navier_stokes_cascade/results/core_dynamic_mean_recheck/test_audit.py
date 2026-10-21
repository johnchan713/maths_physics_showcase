#!/usr/bin/env python3
"""Exact dynamic budget, complete covariance and scientific-scope checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import KEYS,budget,majorant,means,review


class DynamicMeanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=review()
        cls.record=audit.generate()

    def test_independent_current_radial_operator(self):
        self.assertTrue(self.result["checks"]["current_update:full_current_radial_difference"])

    def test_all_current_nonlinear_moment_differences(self):
        for name in ("plain_moment_identity","weighted_moment_identity","theta_current_difference","axial_current_difference"):
            self.assertTrue(self.result["checks"]["current_update:"+name],name)

    def test_old_means_cannot_be_zeroed(self):
        for name in ("old_Haar_means_nonzero","omitted_old_mean_pressure_detected","omitted_old_axial_cross_detected"):
            self.assertTrue(self.result["checks"]["current_update:"+name],name)

    def test_full_covariance_real_wave_rescaling(self):
        for name in ("all_six_covariance_changes_nonzero","complete_wave_rescaling_identity","Drt_changes_target_flux"):
            self.assertTrue(self.result["checks"]["current_update:"+name],name)

    def test_divergence_and_both_moments(self):
        for name in ("old_and_increment_divergence","old_and_new_velocity_moments"):
            self.assertTrue(self.result["checks"]["current_update:"+name],name)

    def test_compact_pressure_and_entire_moving_moment(self):
        for name in ("complete_pressure_both_ends","whole_pressure_axial_moment","moving_compensator_term_nonzero","full_updated_pressure_both_ends","full_updated_pressure_moving_moment"):
            self.assertTrue(self.result["checks"]["current_update:"+name],name)

    def test_coefficient_majorants_satisfy_all_explicit_budgets(self):
        rows={k:v for k,v in self.result["checks"].items() if "norm_bound_" in k or "residual_bound_" in k}
        self.assertEqual(len(rows),16)
        self.assertTrue(all(rows.values()))

    def test_old_mean_and_wave_cross_rate_obstructions(self):
        self.assertTrue(all(self.result["parts"]["correlations"]["checks"].values()))
        rows=self.result["parts"]["correlations"]["rate_rows"]
        self.assertGreater(Q(rows[-1]["old_mean_quadratic_ratio"]),10**10)
        self.assertGreater(Q(rows[-1]["wave_square_ratio"]),10**20)

    def test_zero_update_zero_budget(self):
        n=dict.fromkeys(KEYS,Q(0))
        row=budget(3,Q(1,9),1,2,n)
        self.assertTrue(all(row[k]==0 for k in ("P","Qrem","Jtheta","Jz","Mz","Mtheta")))

    def test_covariance_only_channels(self):
        expected={"Drr":("P","Qrem","Jz"),"Dtt":("P","Qrem","Jz"),
                  "Drz":("P","Qrem","Jz"),"Dtz":("Jtheta",),"Dzz":("Jz",)}
        for key,positive in expected.items():
            n=dict.fromkeys(KEYS,Q(0));n[key]=Q(1)
            row=budget(0,Q(1,5),1,2,n)
            for name in positive:
                self.assertGreater(row[name],0,(key,name))
        n=dict.fromkeys(KEYS,Q(0));n["Drt"]=Q(1)
        row=budget(0,Q(1,5),1,2,n)
        self.assertEqual(row["Drt_kept_for_full_target"],1)
        self.assertTrue(all(row[k]==0 for k in ("P","Qrem","Jtheta","Jz")))

    def test_required_full_ledger_and_rational_inputs(self):
        n=dict.fromkeys(KEYS,Q(0))
        short=dict(n);short.pop("Drt")
        with self.assertRaises(ValueError):
            budget(0,Q(1,5),1,2,short)
        for args in ((-1,Q(1,5),1,2,n),(0,Q(0),1,2,n),(0,Q(1,5),2,1,n)):
            with self.assertRaises(ValueError):
                budget(*args)
        n["B"]=.1
        with self.assertRaises(TypeError):
            budget(0,Q(1,5),1,2,n)

    def test_integer_endpoints_remain_exact(self):
        n=dict.fromkeys(KEYS,Q(1))
        row=budget(2,Q(1,5),1,2,n)
        self.assertTrue(all(isinstance(row[k],Q) for k in ("P","Qrem","Jtheta","Jz","Mz")))

    def test_budget_monotone_in_old_mean_norms(self):
        n=dict.fromkeys(KEYS,Q(1))
        small=budget(2,Q(1,5),1,2,n)
        n["old_v"]=10;n["old_g"]=20
        large=budget(2,Q(1,5),1,2,n)
        for key in ("P","Qrem","Jtheta","Jz"):
            self.assertGreater(large[key],small[key])

    def test_normalized_jet_product_and_derivative_majorants(self):
        Z,T=means.Z,means.T
        f=1+Z+Z**3*T;g=2+T**2+Z*T
        for m in (0,1,3):
            self.assertLessEqual(majorant(f*g,m),majorant(f,m)*majorant(g,m))
            self.assertLessEqual(majorant(f.diff(1),m),(m+1)*majorant(f,m+1))
            self.assertLessEqual(majorant(f.diff(1).diff(1),m),(m+1)*(m+2)*majorant(f,m+2))

    def test_extra_derivative_order_is_explicit(self):
        n=dict.fromkeys(KEYS,Q(1))
        row=budget(4,Q(1,5),1,2,n)
        self.assertEqual(row["input_order"],6)
        self.assertEqual(row["residual_order"],3)
        self.assertEqual(budget(0,Q(1,5),1,2,n)["Mz"],None)

    def test_both_pressure_tails_retained(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))

    def test_historical_source_bindings(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)

    def test_no_uniform_scale_or_full_residual_promotion(self):
        for key in ("generic_correlations_are_new_square_order","norm_gate_verified_uniformly_on_all_actual_scales","actual_five_map_or_conditional_cycle_refuted_by_control","pointwise_full_residual_improvement_or_blowup_proved","polynomial_diagnostics_are_continuum_suprema"):
            self.assertFalse(self.record["flags"][key])

    def test_scope_tampering_rejected(self):
        record=deepcopy(self.record)
        record["flags"]["norm_gate_verified_uniformly_on_all_actual_scales"]=True
        with self.assertRaises(ValueError):
            audit.validate(record)

    def test_frozen_record_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_output_overwrite_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"evidence.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):
                audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":
    unittest.main()
