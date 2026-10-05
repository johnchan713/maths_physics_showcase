#!/usr/bin/env python3
"""Meaningful controls for width losses, angular pressure and unsupported promotion."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest
from bounds import AMIN,HERE,radius,repair_budget,scalar_checks,ledger
from review import bump_review,moment_pressure_review
from audit import STATUS,FLAGS,RESOLVED,REMAINING,validate,provenance,expected_protocol


class RepairTests(unittest.TestCase):
    def test_positive_coefficient_neighborhood_uses_exact_small_numbers(self):
        b=repair_budget(AMIN,radius(AMIN))
        self.assertGreater(b['coefficient_norm'],0)
        self.assertLess(b['relative_swirl_error'],F(1,2))
        self.assertLess(b['field_C1'],1)
        self.assertLess(b['shear_error'],AMIN**6*b['coefficient_norm'])

    def test_comparison_neighborhood_rejects_oversized_coefficients(self):
        with self.assertRaises(ValueError): repair_budget(AMIN,2*radius(AMIN))
        with self.assertRaises(ValueError): repair_budget(AMIN,-radius(AMIN))
        with self.assertRaises(ValueError): radius(100)

    def test_bump_width_and_angular_denominator_controls(self):
        r=bump_review()
        self.assertTrue(all(r['checks'].values()),r['checks'])

    def test_five_moments_original_source_ODEs_and_omissions(self):
        r=moment_pressure_review()
        self.assertTrue(all(r['checks'].values()),r['checks'])

    def test_outward_constants_and_all_frequency_thresholds(self):
        for digits in (80,110):
            r=scalar_checks(digits)
            self.assertTrue(all(r['checks'].values()),r['checks'])
            self.assertFalse(r['huge_frequency_materialized'])

    def test_state_bound_does_not_invent_a_pressure_derivative(self):
        r=ledger()
        self.assertEqual(r['pressure_state_derivative_order'],0)
        self.assertEqual(r['moment_angular_order'],1)
        self.assertEqual(r['state_coefficient'],327680)
        self.assertTrue(all(r['checks'].values()))

    def test_frozen_dependencies_and_protocol(self):
        self.assertTrue(all(provenance().values()))
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),expected_protocol())


class ScopeTests(unittest.TestCase):
    def record(self):
        return deepcopy(dict(status=STATUS,flags=FLAGS,resolved_obligations=RESOLVED,
                             remaining_obligations=REMAINING,gates={'valid':True}))

    def test_frequency_acceptance_remains_conditional(self):
        validate(self.record())
        self.assertTrue(FLAGS['relies_on_inherited_analytic_construction'])
        self.assertFalse(FLAGS['independent_foundation_review_completed'])

    def test_physical_and_blowup_promotions_rejected(self):
        for key in ('full_physical_stress_wave_construction_verified',
                    'full_PDE_corrections_verified','smooth_force_verified',
                    'blowup_verified','independent_foundation_review_completed'):
            r=self.record();r['flags'][key]=True
            with self.assertRaises(ValueError): validate(r)

    def test_remaining_obligations_and_failed_checks_cannot_be_hidden(self):
        r=self.record();r['remaining_obligations']=[]
        with self.assertRaises(ValueError): validate(r)
        r=self.record();r['gates']['valid']=False
        with self.assertRaises(ValueError): validate(r)


if __name__=='__main__': unittest.main()
