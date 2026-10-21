#!/usr/bin/env python3
"""Exact implicit-chart, moving support and compensator scope controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import Jet,c_norm,dynamic,geometry,implicit_jet,means,reserved_budget,review


class ReservedGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def test_complex_contraction_constants(self):
        self.assertLess(F(4,125),1)
        self.assertEqual(F(49,256)+F(1,125),F(6381,32000))
        self.assertLess(F(6381,32000),F(1,4))

    def test_defining_equation_and_Cauchy_bounds_all_fixture_centers(self):
        checks=self.result["parts"]["coordinate"]["checks"]
        self.assertTrue(all(checks.values()))

    def test_exact_first_derivatives_and_polar_points(self):
        for eta in (F(-1),F(0),F(1)):
            h=F(1,200);t=implicit_jet(h,eta,3);L=1-2*h*eta**2
            self.assertEqual(t.c.get((1,0),0),2*eta/L)
            self.assertEqual(t.c.get((0,1),0),1/L)

    def test_omitted_implicit_divisor_is_detected(self):
        t=implicit_jet(F(1,200),F(1),2)
        self.assertNotEqual(t.c[(0,1)],1)
        self.assertNotEqual(t.c[(1,0)],2)

    def test_tiny_fixed_h_uses_only_rational_coordinate_jets(self):
        t=implicit_jet(F(1,2**80),F(1),3)
        self.assertEqual(t.c[(0,1)],1/(1-F(2,2**80)))
        self.assertTrue(all(isinstance(v,F) for v in t.c.values()))

    def test_invalid_chart_data_rejected(self):
        for h,eta in ((F(0),F(1)),(F(1,2),F(1)),(F(1,200),F(2))):
            with self.assertRaises(ValueError):implicit_jet(h,eta,2)
        with self.assertRaises(TypeError):implicit_jet(.005,F(1),2)
        with self.assertRaises(ValueError):Jet(-1)

    def test_series_power_requires_zero_constant_and_same_order(self):
        with self.assertRaises(ValueError):Jet(2,{(0,0):1}).power_one(F(1,2))
        with self.assertRaises(ValueError):_ = Jet(2)+Jet(3)
        with self.assertRaises(ValueError):_ = Jet(2)*Jet(3)

    def test_compensator_integrals_on_moving_intervals(self):
        self.assertTrue(all(self.result["parts"]["radial_compensator"]["checks"].values()))

    def test_compensator_norm_bound_for_all_finite_orders(self):
        self.assertEqual(c_norm(0,F(3,7)),F(12,7))
        self.assertEqual(c_norm(1,F(3,7)),F(300,7))
        for n in (1,3,10):self.assertGreater(c_norm(n,F(3,7)),c_norm(n-1,F(3,7)))
        with self.assertRaises(ValueError):c_norm(2,F(-1))

    def test_signed_unit_mass_does_not_bound_absolute_second_moment(self):
        R=means.R;b=(R-1)**4*(2-R)**4;hat=b/b.integral().scalar()
        signed=hat+100*(R-F(3,2))*hat
        self.assertEqual(signed.integral(),1)
        self.assertGreater((R**2*signed).integral().scalar()/2,F(2))

    def test_logarithmic_gap_certificate_keeps_actual_width(self):
        g=geometry(2,3,F(1,2),16,F(1,2**20))
        self.assertEqual(g["left_log_gap_lower"],F(6,5))
        self.assertEqual(g["right_log_gap_lower"],F(46,41))
        other=geometry(2,3,F(1,2),16,F(1,2**10))
        self.assertGreater(other["log_weight_comparison"],g["log_weight_comparison"])

    def test_hull_must_stay_strictly_inside_active_annulus(self):
        for args in ((2,3,2,16,1),(2,3,F(1,2),F(9,2),1),(3,2,F(1,2),16,1)):
            with self.assertRaises(ValueError):geometry(*args)

    def test_moving_support_and_rectangular_hull_distinguished(self):
        self.assertTrue(all(self.result["parts"]["geometry_budget"]["checks"].values()))
        self.assertEqual(F(self.result["parts"]["geometry_budget"]["outside_rectangular_point_X"]),72)

    def test_dynamic_geometry_and_compensator_substitution(self):
        n={k:F(1) for k in dynamic.KEYS if k!="c_rho"}
        for k in ("Drr","Drt","Drz","Dtt","Dtz","Dzz"):n[k]=F(0)
        got=reserved_budget(2,F(1,7),2,3,F(3,7),n)
        direct=dict(n);direct["c_rho"]=c_norm(4,F(3,7))
        self.assertEqual(got,dynamic.budget(2,F(1,7),1,6,direct))

    def test_reserved_budget_refuses_changing_wave_or_incomplete_ledger(self):
        n={k:F(0) for k in dynamic.KEYS if k!="c_rho"}
        n["Drr"]=F(1)
        with self.assertRaises(ValueError):reserved_budget(2,F(1,7),2,3,F(3,7),n)
        n.pop("Drt")
        with self.assertRaises(ValueError):reserved_budget(2,F(1,7),2,3,F(3,7),n)

    def test_both_full_background_pressure_tails_retained(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))

    def test_declared_historical_bindings_pass(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)

    def test_coordinate_result_does_not_promote_heat_or_all_field_bounds(self):
        for key in ("rectangular_union_hull_has_positive_weight_floor","coordinate_analytic_extension_is_heat_profile_analytic_extension","all_actual_generated_field_norms_certified","reserved_hull_gate_covers_wider_changing_wave_supports","complete_residual_improvement_or_blowup_proved","order_six_diagnostics_prove_all_orders"):
            self.assertFalse(self.record["flags"][key])

    def test_scope_tampering_rejected(self):
        record=deepcopy(self.record);record["flags"]["all_actual_generated_field_norms_certified"]=True
        with self.assertRaises(ValueError):audit.validate(record)

    def test_frozen_evidence_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_output_overwrite_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"evidence.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
