#!/usr/bin/env python3
"""Exact controls for coefficient supports and reserved pressure."""

from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import (
    axial_cut, coefficient_bound, conservative, cutoff_flux, denormalize,
    normalize, pressure_departure, pressure_floor, radial_flux,
    reciprocal_norm, review, support_edges, time_cut,
)


class CoefficientSupportRecheckTests(unittest.TestCase):
    def test_ordinary_rows_roundtrip_with_nonconstant_amplitude(self):
        raw = tuple(map(Q, (2, 3, 5, 7, 11)))
        for eta in (Q(-1), Q(0), Q(2, 5), Q(1)):
            scales = (Q(2), Q(3, 5), Q(7) / (1 + eta * eta), Q(1, 1000))
            self.assertEqual(denormalize(normalize(raw, *scales), *scales), raw)
        self.assertNotEqual(normalize(raw, 2, Q(3, 5), 7, Q(1, 1000)), normalize(raw, 2, Q(3, 5), Q(7, 2), Q(1, 1000)))

    def test_zero_lambda_rejected(self):
        with self.assertRaises(ValueError):
            normalize((1,) * 5, 1, 1, 1, 0)

    def test_current_order_squares_have_later_grade(self):
        for n in (1, 2, 3, 19, 100):
            pairs = [(i, n - i) for i in range(n + 1)]
            self.assertEqual([p for p in pairs if n in p], [(0, n), (n, 0)])
            self.assertNotIn((n, n), pairs)

    def test_full_reciprocal_jet_cost(self):
        self.assertEqual([reciprocal_norm(k, 7) for k in (0, 1, 2, 1000000)], [Q(2, 7), Q(4, 7), Q(5, 7), Q(5, 7)])
        with self.assertRaises(ValueError):
            reciprocal_norm(-1, 7)

    def test_coefficient_bound_keeps_lambda_loss(self):
        args = ((2, 3, 5, 7, 11), 8, 2, Q(3, 5), 7)
        self.assertEqual(coefficient_bound(*args, Q(1, 1000)), 3200000)
        self.assertEqual(coefficient_bound(*args, Q(1, 2000)), 6400000)

    def test_bump_support_before_reserved_interval(self):
        centers = tuple(map(Q, ("1/2", "1", "3/2", "2", "5/2")))
        self.assertLess(max(c + Q(1, 20) for c in centers), 3)
        self.assertLess(3, 6)

    def test_leading_accumulated_flux_remains(self):
        for eta in map(Q, ("-1", "0", "1/3", "1")):
            self.assertEqual(radial_flux(0, 3 * eta, Q(3), Q(10), eta, Q(1, 200), 0), -3)

    def test_all_positive_velocity_and_potential_cutoff_fluxes_vanish(self):
        for n in (1, 2, 3, 7, 1000000):
            V = radial_flux(0, 0, 0, Q(10), Q(2, 5), Q(1, 200), n)
            self.assertEqual(cutoff_flux(Q(5, 32), Q(-27, 32), V, 0, Q(2, 5), Q(1, 200)), 0)

    def test_M_jet_cannot_be_replaced_by_value_match(self):
        self.assertNotEqual(radial_flux(0, 0, Q(1), Q(10), Q(0), Q(1, 200), 1), 0)

    def test_omitting_cutoff_curl_changes_nonzero_memory(self):
        self.assertNotEqual(cutoff_flux(0, Q(1), 0, Q(1), Q(1, 2), Q(1, 200)), 0)

    def test_pressure_and_stress_support_exceptions(self):
        self.assertEqual(support_edges(1), {"velocity": "Xvel", "pressure": "Xv", "stress": "Xb"})
        self.assertEqual(support_edges(2), {"velocity": "Xvel", "pressure": "Xv", "stress": "Xv"})
        self.assertEqual(set(support_edges(1000000).values()), {"Xvel"})

    def test_pressure_floor_survives_arbitrary_constants(self):
        a, b, mu, chi, scale = Q(10), Q(30), Q(3), Q(3, 7), Q(1, 5)
        floor = pressure_floor(a, b, mu, chi, scale)
        for K1, K2 in ((Q(7), Q(-13)), (Q(-100), Q(300)), (Q(0), Q(0))):
            p0 = pressure_departure(a, mu, K1, K2, chi, scale, Q(2, 3))
            p1 = pressure_departure(b, mu, K1, K2, chi, scale, Q(2, 3))
            self.assertEqual(p1 - p0, 2 * floor)
            self.assertGreaterEqual(max(abs(p0), abs(p1)), floor)

    def test_midpoint_constant_saturates_pressure_floor(self):
        a, b, mu = Q(10), Q(30), Q(3)
        mid = mu * mu * (1 / a + 1 / b) / 8
        floor = pressure_floor(a, b, mu, Q(1), Q(1))
        self.assertEqual([pressure_departure(x, mu, mid, 0, Q(1), Q(1), 0) for x in (a, b)], [-floor, floor])

    def test_inactive_first_cutoff_has_zero_obstruction(self):
        self.assertEqual(pressure_floor(Q(10), Q(30), Q(3), Q(0), Q(1, 5)), 0)
        with self.assertRaises(ValueError):
            pressure_floor(Q(30), Q(10), Q(3), Q(1), Q(1))

    def test_inverse_radius_piece_is_stationary_when_uncut(self):
        X, mu, eta, h = Q(10), Q(3), Q(2, 5), Q(1, 200)
        P, PX = -mu * mu / (4 * X), mu * mu / (4 * X * X)
        self.assertEqual(axial_cut(-1, P, PX, 0, X, eta, h, 1, 0), 0)
        self.assertEqual(time_cut(-1, P, PX, 0, X, eta, h, 1, 0), 0)
        for q in (Q(1, 7), Q(1), Q(13)):
            r_squared = 2 * q * X
            self.assertEqual(P / q, -mu * mu / (2 * r_squared))

    def test_cutoff_transition_matches_direct_physical_time_derivative(self):
        # chi(cq)=(cq)^2: differentiate -a*c^2*q^2/(q*X)
        # at fixed physical r, using q_t=-1/L and r^2=2qX.
        q, c, X, eta, h, a = Q(2), Q(3), Q(10), Q(2, 5), Q(1, 200), Q(9, 4)
        L = 1 - 2 * h * eta * eta
        chi = (c * q) ** 2
        actual = time_cut(-1, -a / X, a / (X * X), 0, X, eta, h, chi, 2 * chi) / (q * q)
        self.assertEqual(actual, 2 * a * c * c / (X * L))
        self.assertNotEqual(actual, 0)

    def test_cutoff_transition_matches_direct_physical_axial_derivative(self):
        q, c, X, eta, h, a = Q(2), Q(3), Q(10), Q(2, 5), Q(1, 200), Q(9, 4)
        L = 1 - 2 * h * eta * eta
        chi = (c * q) ** 2
        # Factor q^-D from both physical expressions to keep exact rationals.
        actual = axial_cut(-1, -a / X, a / (X * X), 0, X, eta, h, chi, 2 * chi) / q
        self.assertEqual(actual, -4 * a * c * c * eta * q / (X * L))
        self.assertNotEqual(actual, 0)

    def test_angular_pressure_constant_derivative_is_not_dropped(self):
        self.assertNotEqual(axial_cut(-1, 0, 0, Q(1), Q(10), 0, Q(1, 200), 1, 0), 0)

    def test_nonlocal_second_pressure_control(self):
        self.assertEqual(conservative.plateau_review()["manufactured_second_pressure_plateau"], "153/10")

    def test_radially_constant_second_pressure_has_physical_derivatives(self):
        h, eta, K2 = Q(1, 200), Q(2, 5), Q(7)
        b = -1 + 2 * h
        self.assertNotEqual(time_cut(b, K2, 0, 0, Q(10), eta, h, 1, 0), 0)
        self.assertNotEqual(axial_cut(b, K2, 0, 0, Q(10), eta, h, 1, 0), 0)

    def test_nonzero_lower_viscosity_in_order_seven_identity(self):
        direct, expected, viscous = conservative.assemble(7, Q(1, 200), Q(1, 3), Q(2))
        self.assertEqual(direct, expected)
        self.assertTrue(any(viscous))

    def test_all_review_gates_and_dependency_sources(self):
        record = audit.generate()
        self.assertTrue(all(review()["checks"].values()))
        self.assertTrue(all(record["gates"].values()))
        self.assertGreater(len(record["inherited_source_checks"]), 50)

    def test_false_scope_and_erased_frontier_are_rejected(self):
        record = audit.generate()
        for key, value in audit.FLAGS.items():
            if not value:
                changed = deepcopy(record)
                changed["flags"][key] = True
                with self.assertRaises(ValueError):
                    audit.validate(changed)
        changed = deepcopy(record)
        changed["remaining_obligations"].pop()
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_frozen_record_and_fresh_replay(self):
        audit.validate(json.loads((audit.HERE / "evidence.json").read_text()))
        with tempfile.TemporaryDirectory(dir=audit.PROJECT) as directory:
            path = Path(directory) / "fresh.json"
            audit.main(["--verify-record", str(audit.HERE / "evidence.json"), "--output", str(path)])
            audit.validate(json.loads(path.read_text()))
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])


if __name__ == "__main__":
    unittest.main()
