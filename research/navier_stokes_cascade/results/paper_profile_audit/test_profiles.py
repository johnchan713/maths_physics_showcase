#!/usr/bin/env python3
"""Independent calculations and deliberate failure controls; no proof claims."""
import copy
from fractions import Fraction
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from profiles import (Jet, ManufacturedProfile, Poly, coordinates, core_coefficients,
                      core_comparison, cylindrical_and_stress, field,
                      finite_difference_residual, full_residual, heat_derivatives,
                      heat_exterior, heat_pressure, heat_reference, normalized_error,
                      physical_point, q_value)
from run_audit import HERE, STATUS, audit, compare_reproduction, main, validate, verify_record

PROTOCOL = json.loads((HERE/"protocol.json").read_text())


class DerivativeTests(unittest.TestCase):
    def test_polynomial_jet(self):
        x, y = Jet.variable(2., 0), Jet.variable(-3., 1)
        result = x*x*y + 3*y*y
        self.assertEqual(result.v, 15.)
        np.testing.assert_array_equal(result.g, [-12., -14., 0., 0.])
        np.testing.assert_array_equal(result.H[:2, :2], [[-6., 4.], [4., 6.]])

    def test_quotient_jet(self):
        x, y = Jet.variable(2., 0), Jet.variable(4., 1)
        result = x/y
        self.assertEqual(result.v, .5)
        np.testing.assert_array_equal(result.g, [.25, -.125, 0., 0.])
        np.testing.assert_array_equal(result.H[:2, :2], [[0., -.0625], [-.0625, .0625]])

    def test_axis_polynomial_power(self):
        x = Jet.variable(0., 0)
        self.assertEqual((x**0).v, 1.)
        self.assertEqual((x**2).H[0, 0], 2.)
        self.assertEqual((x**3).H[0, 0], 0.)
        for power in (.5, -1):
            with self.assertRaises(ValueError):
                x**power

    def test_polynomial_integral_and_regular_average(self):
        p = Poly({(0, 0): 2., (1, 1): 3., (2, 0): 6.})
        self.assertEqual(p.integral().derivative(0).c, p.c)
        self.assertEqual(p.average()(0., .7), 2.)
        for X in (.2, 1.):
            self.assertAlmostEqual(p.integral()(X, .3)/X, p.average()(X, .3))

    def test_q_root_and_z_symmetry(self):
        for h in (.005, .2):
            for z in (0., .3, -1.):
                q = q_value(z, .7, h)
                self.assertAlmostEqual(q-z*z*q**(2*h), .3, places=13)
                self.assertEqual(q, q_value(-z, .7, h))

    def test_q_implicit_hessian_by_differencing(self):
        point, h, step = np.array([.2, -.3, .4, .6]), .005, 1e-4
        q = coordinates(point, h)[4]
        for i in (2, 3):
            plus, minus = point.copy(), point.copy()
            plus[i] += step
            minus[i] -= step
            qp, qm = coordinates(plus, h)[4], coordinates(minus, h)[4]
            self.assertLess(abs((qp.v-qm.v)/(2*step)-q.g[i]), 2e-9)
            np.testing.assert_allclose((qp.g-qm.g)/(2*step), q.H[:, i], rtol=1e-7, atol=2e-9)
        np.testing.assert_array_equal(q.H, q.H.T)

    def test_coordinate_domains(self):
        for point in ([0, 0, 0, 1], [0, 0, 0, -.1], [math.nan, 0, 0, .5]):
            with self.assertRaises(ValueError):
                coordinates(point, .005)
        for h in (0, .5, math.nan):
            with self.assertRaises(ValueError):
                q_value(0., .2, h)
        for point in ((-1, 0, 1, 0), (.2, 1, 1, 0), (.2, 0, 1, math.inf)):
            with self.assertRaises(ValueError):
                physical_point(*point, .005)

    def test_normalized_error_refuses_invalid(self):
        for a, b in (([1], [math.nan]), ([], []), ([1, 2], [1])):
            with self.assertRaises(ValueError):
                normalized_error(a, b)
        with self.assertRaises(ValueError):
            normalized_error(1., 1., floor=0)


class EquationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.profile = ManufacturedProfile(PROTOCOL["fixtures"])
        cls.h = PROTOCOL["fixtures"]["h"]
        cls.point = physical_point(.2, -.6, .7, .4, cls.h)

    def test_manufactured_profile_is_not_solution(self):
        u, p, _ = field(self.point, self.h, self.profile)
        residual, scale, _ = full_residual(u, p)
        self.assertGreater(np.max(np.abs(residual))/scale, .01)
        self.assertLess(abs(sum(u[i].g[i] for i in range(3))), 1e-12)

    def test_cartesian_cylindrical_and_stress(self):
        u, p, _ = field(self.point, self.h, self.profile)
        residual, _, _ = full_residual(u, p)
        c = cylindrical_and_stress(self.point, self.h, self.profile)
        self.assertLess(normalized_error(residual, c["cartesian"]), 1e-12)
        self.assertLess(normalized_error(c["leading_tangential"], c["stress_rhs"]), 1e-12)
        self.assertGreater(normalized_error(c["cylindrical"][1:], c["stress_rhs"]), .01)

    def test_finite_difference_convergence(self):
        u, p, _ = field(self.point, self.h, self.profile)
        residual, scale, _ = full_residual(u, p)
        errors = [np.max(np.abs(finite_difference_residual(self.point, self.h, self.profile, step)-residual))/scale
                  for step in (.02, .01)]
        self.assertGreater(errors[0]/errors[1], 10.)
        self.assertLess(errors[1], 1e-6)

    def test_axis_regular_extension(self):
        u, p, _ = field([0., 0., .3, .5], self.h, self.profile)
        self.assertEqual(u[0].v, 0.)
        self.assertEqual(u[1].v, 0.)
        self.assertLess(abs(sum(u[i].g[i] for i in range(3))), 1e-12)
        self.assertTrue(np.isfinite(full_residual(u, p)[0]).all())

    def test_comparison_exact_recurrence(self):
        coefficients = core_coefficients()
        self.assertEqual(coefficients[:4], [Fraction(1), Fraction(-1, 4), Fraction(1, 48), Fraction(-1, 1152)])
        for n in range(len(coefficients)-1):
            self.assertEqual(2*(n+1)*(n+2)*coefficients[n+1]+coefficients[n], 0)
        # A finite truncation is not an exact solution of the full ODE.
        self.assertNotEqual(coefficients[-1], 0)

    def test_comparison_lower_bound(self):
        t = Fraction(41, 20)
        lower = 1-t/2+t*t/12-t**3/144
        self.assertEqual(lower, Fraction(305719, 1152000))
        for numerator in range(42):
            x = Fraction(numerator, 20)
            self.assertEqual(-Fraction(1, 2)+x/6-x*x/48, -((x-4)**2+8)/48)
        for z in PROTOCOL["core_comparison"]["z"]:
            value, reference, bound = core_comparison(z)
            self.assertGreater(value, float(lower))
            self.assertLessEqual(value, 1.)
            self.assertLess(abs(value-reference), 1e-12)
            self.assertLess(bound, 1e-40)

    def test_core_invalid_inputs(self):
        for z in (-.1, 4.2, math.nan):
            with self.assertRaises(ValueError):
                core_comparison(z)
        for terms in (3, 101, 5.5):
            with self.assertRaises(ValueError):
                core_coefficients(terms)

    def test_heat_endpoint_derivatives(self):
        values, _ = heat_derivatives(.005, 0.)
        np.testing.assert_allclose(values, [heat_reference(.005, 0., j) for j in range(4)], rtol=1e-11, atol=1e-13)
        self.assertAlmostEqual(values[0], 1., places=12)

    def test_heat_full_equation_and_negative_control(self):
        result = heat_exterior([.7, .2, 0., .2], .005)
        self.assertLess(np.max(np.abs(result["residual"]))/result["scale"], 1e-11)
        self.assertGreater(abs(result["without_swirl_curvature"])/result["scale"], 1e-4)
        self.assertAlmostEqual(result["without_swirl_curvature"], result["expected_curvature_error"], places=11)

    def test_heat_pressure_gradient_from_separate_integral(self):
        radius, time, h, step = .9, .6, .005, .0002
        values = [heat_pressure(radius+j*step, time, h) for j in (-2, -1, 1, 2)]
        derivative = (values[0]-8*values[1]+8*values[2]-values[3])/(12*step)
        direct = heat_exterior([radius, 0., 0., time], h)["pressure_gradient"][0]
        self.assertLess(normalized_error(derivative, direct), 1e-9)

    def test_heat_domains_and_final_time_away_axis(self):
        for point in ([0, 0, 0, .5], [1, 0, 0, 1.1], [1, 0, 0, math.nan]):
            with self.assertRaises(ValueError):
                heat_exterior(point, .005)
        for Z in (-1, math.inf):
            with self.assertRaises(ValueError):
                heat_derivatives(.005, Z)
        result = heat_exterior([1., 0., .3, 1.], .005)
        self.assertLess(np.max(np.abs(result["residual"]))/result["scale"], 1e-11)


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = audit(PROTOCOL)

    def test_complete_audit(self):
        validate(self.report, PROTOCOL)
        self.assertEqual(self.report["status"], STATUS)
        self.assertEqual(len(self.report["manufactured_samples"]), 54)
        self.assertEqual(self.report["axis_samples"], 9)
        self.assertEqual(len(self.report["heat_exterior"]), 6)

    def test_missing_check_is_rejected(self):
        report = copy.deepcopy(self.report)
        del report["checks"]["heat-full-residual"]
        with self.assertRaises(ValueError):
            validate(report, PROTOCOL)

    def test_missing_or_nonfinite_raw_samples_rejected(self):
        report = copy.deepcopy(self.report)
        report["manufactured_samples"].pop()
        with self.assertRaises(ValueError):
            validate(report, PROTOCOL)
        report = copy.deepcopy(self.report)
        report["manufactured_samples"][0]["full_cartesian_residual"][0] = math.nan
        with self.assertRaises(ValueError):
            validate(report, PROTOCOL)

    def test_changed_physical_evidence_rejected(self):
        report = copy.deepcopy(self.report)
        report["manufactured_samples"][0]["full_cartesian_residual"][0] += 1.
        with self.assertRaisesRegex(ValueError, "not reproduced"):
            compare_reproduction(report, self.report, PROTOCOL)
        compare_reproduction(self.report, self.report, PROTOCOL)

    def test_nonfinite_and_negative_error_rejected(self):
        for bad in (math.nan, math.inf, -1., True):
            report = copy.deepcopy(self.report)
            report["checks"]["heat-full-residual"]["value"] = bad
            with self.assertRaises(ValueError):
                validate(report, PROTOCOL)

    def test_changed_limit_rejected(self):
        report = copy.deepcopy(self.report)
        report["checks"]["heat-full-residual"]["limit"] = 100.
        with self.assertRaises(ValueError):
            validate(report, PROTOCOL)

    def test_failed_negative_controls_cannot_pass(self):
        for name in PROTOCOL["required_negative_controls"]:
            report = copy.deepcopy(self.report)
            report["checks"]["detect-"+name]["value"] = 0.
            report["status"] = "building-block-audit-failed"
            with self.assertRaisesRegex(ValueError, "Audit gates failed"):
                validate(report, PROTOCOL)

    def test_promotion_and_obligation_erasure_rejected(self):
        for key, value in (("full_proof_verified", True), ("new_candidate", True),
                           ("status", "proof-verified"), ("unverified_obligations", [])):
            report = copy.deepcopy(self.report)
            report[key] = value
            with self.assertRaises(ValueError):
                validate(report, PROTOCOL)

    def test_source_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"evidence.json"
            path.write_text(json.dumps({"source_sha256": {"profiles.py": "wrong"}}))
            with self.assertRaisesRegex(ValueError, "source files"):
                verify_record(path, PROTOCOL, {"profiles.py": "expected"})

    def test_existing_output_and_wrong_pdf_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"untouched.json"
            path.write_text("preserve")
            with self.assertRaises(SystemExit):
                main(["--output", str(path)])
            with self.assertRaisesRegex(ValueError, "PDF bytes"):
                main(["--output", str(Path(directory)/"new.json"), "--paper", str(path)])
            self.assertEqual(path.read_text(), "preserve")
            self.assertFalse((Path(directory)/"new.json").exists())


if __name__ == "__main__":
    unittest.main()
