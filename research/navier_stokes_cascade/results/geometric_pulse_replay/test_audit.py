#!/usr/bin/env python3
"""Radius-aware analytic certificates, independent ODE checks and scope controls."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
import numpy as np
from scipy.integrate import quad
from estimates import (RADIUS, all_later_threshold, covariance_interval,
                       derivative_and_tail_review, derivative_degrees,
                       exact_scale_review, flat_tail_threshold, gaussian_constants,
                       length_interval, minimum_pulse_length, pulse_barrier,
                       radius, radius_power_ledger, real, threshold_review)
from diagnostics import (Fixture, beta, covariance_review, failed_shortcuts_review,
                         frame_review, harmonic_factor_review, lifted_source_review,
                         psi, sourced_review, transverse_cutoff)
from audit import (FLAGS, HERE, compare_records, expected_protocol, generate,
                   inherited_provenance, main, refinement_review, validate)


class LengthAndDomainTests(unittest.TestCase):
    def test_radius_is_the_exact_constructive_global_value(self):
        self.assertEqual(RADIUS,F(1,6057072))
        self.assertTrue(all(exact_scale_review()['checks'].values()))

    def test_pulse_length_interval_retains_radius(self):
        with mp.workdps(60):
            lo,hi=length_interval(51600)
            self.assertEqual(lo,2*real(RADIUS)*51600**2)
            self.assertLess(hi,51600**2)
            self.assertGreater(hi,lo)

    def test_invalid_radius_band_and_threshold_are_rejected(self):
        for r in (0,-1,F(1,25),float('inf')):
            with self.assertRaises(ValueError):
                radius(r)
        for ell in (15,16.):
            with self.assertRaises(ValueError):
                length_interval(ell)
        for h,r,length,first in ((.02,RADIUS,25,16),(.005,RADIUS,0,16),(.005,RADIUS,25,0)):
            with self.assertRaises(ValueError):
                all_later_threshold(h,r,length,first)

    def test_geometry_validity_does_not_accept_short_analytic_pulses(self):
        self.assertLess(length_interval(16)[1],1)
        with self.assertRaises(ValueError):
            pulse_barrier(1,1,0,length_interval(16)[1])

    def test_sufficient_threshold_applies_to_every_later_band(self):
        row=threshold_review()
        self.assertTrue(all(row['checks'].values()))
        for value in row['rows']:
            lo,_=length_interval(2*value['first_band'],value['radius'])
            self.assertGreaterEqual(lo,4*value['minimum_length'])

    def test_smaller_radius_increases_length_threshold_at_fixed_h(self):
        a=all_later_threshold(F(1,200),RADIUS,1000)
        b=all_later_threshold(F(1,200),RADIUS/4,1000)
        self.assertGreater(b['pulse_length_band'],a['pulse_length_band'])
        self.assertEqual(a['exponential_band'],b['exponential_band'])

    def test_cone_threshold_and_riccati_length_are_both_retained(self):
        length=minimum_pulse_length('.4',2,3,2,1,'.2')
        self.assertGreaterEqual(length,40)
        self.assertLessEqual(3/mp.sqrt(length)*(mp.mpf('.5')+1),mp.mpf('.2')/mp.mpf('4.2')+mp.mpf('1e-15'))


class PulseAndCovarianceTests(unittest.TestCase):
    def test_actual_phase_formula_at_admissible_radius_has_length_scaled_error(self):
        self.assertTrue(all(frame_review()['checks'].values()))

    def test_riccati_barrier_is_strict_and_radial_lower_bound_positive(self):
        for length in (40,100,1000):
            row=pulse_barrier('.4',2,'.5',length)
            self.assertLess(row['outward_derivative_upper'],0)
            self.assertGreater(row['radial_lower'],0)
            self.assertLessEqual(row['ratio_bound'],mp.mpf('.125'))

    def test_pulse_bounds_do_not_hide_a_radius_factor_in_the_forward_growth(self):
        a=pulse_barrier('.4',2,'.5',100)
        b=pulse_barrier('.4',2,'.5',1000)
        self.assertEqual(a['log_forward_propagator_factor'],b['log_forward_propagator_factor'])
        self.assertLess(a['log_amplitude_bound'],3)

    def test_Gaussian_constants_bound_the_reference_second_log_derivative(self):
        with mp.workdps(50):
            g=gaussian_constants('.3','.5','1.5')
            for lam in (mp.mpf('.3'),mp.mpf('.5')):
                for y in (mp.mpf('.75'),mp.mpf('1.5'),mp.mpf('2.25')):
                    value=mp.mpf('1.5')*lam*y*((1+y*y)**(-mp.mpf('1.5'))+2*(1+mp.mpf('1.5')**2)**(-mp.mpf('1.5')))
                    self.assertLessEqual(g['lower'],value)
                    self.assertGreaterEqual(g['upper'],value)

    def test_higher_harmonic_scalar_factor_and_negative_harmonic(self):
        self.assertTrue(all(harmonic_factor_review()['checks'].values()))

    def test_full_Haar_covariance_is_positive_at_the_global_radius(self):
        row=covariance_review()
        self.assertTrue(all(row['checks'].values()))
        for r in row['rows']:
            self.assertLess(r['h'],1e-13)
            self.assertGreater(r['h'],0)

    def test_normalized_covariance_scales_with_radius_squared_and_length(self):
        for row in covariance_review()['rows']:
            self.assertTrue(.3<row['h_times_sqrt_length_over_radius_squared']<.6)
            self.assertTrue(row['full_pulse_interval_integrated'])

    def test_exact_radius_power_ledger_retains_amplitude_compensation(self):
        p=radius_power_ledger()
        self.assertEqual(p['haar_covariance_r'],F(3,2))
        self.assertEqual(p['covariance_inverse_r'],F(-3,2))
        self.assertEqual(p['primary_root_r'],F(-3,4))
        self.assertEqual(p['squared_primary_root_times_haar_r'],0)

    def test_continuum_covariance_interval_has_exact_radius_squared_scaling(self):
        a=covariance_interval(RADIUS,100,.1,2,.5,2)
        b=covariance_interval(RADIUS/2,100,.1,2,.5,2)
        self.assertAlmostEqual(float(a['lower']/b['lower']),4)
        self.assertAlmostEqual(float(a['upper']/b['upper']),4)
        with self.assertRaises(ValueError):
            covariance_interval(RADIUS,24,.1,2,.5,2)

    def test_transverse_cutoff_area_is_not_an_implicit_unit(self):
        integral,_=quad(lambda x:beta(x)**2,-1,1)
        self.assertGreater(integral,float(mp.exp(-mp.mpf(2)/3)))
        self.assertLess(integral,2)
        self.assertLess(float(RADIUS)/2*integral,1e-6)

    def test_transverse_cutoff_is_compactly_inside_the_original_open_rectangle(self):
        r=float(RADIUS)
        self.assertEqual(transverse_cutoff(0,r),1)
        for factor in (-1,-.75,-.5,.5,.75,1):
            self.assertEqual(transverse_cutoff(factor*r,r),0)

    def test_tighter_ODE_and_full_pulse_quadrature_reproduce_covariance(self):
        row=refinement_review(covariance_review(),covariance_review(True))
        self.assertTrue(all(row['checks'].values()))


class SourcedInverseTests(unittest.TestCase):
    def test_independent_two_frame_and_three_component_solutions_agree(self):
        row=sourced_review()
        self.assertTrue(row['checks']['independent_projected_and_frame_solutions_agree'])

    def test_moving_normal_constraint_and_complete_pressure_balance(self):
        row=sourced_review()
        self.assertTrue(row['checks']['moving_normal_constraint_preserved'])
        self.assertTrue(row['checks']['full_principal_pressure_and_viscosity_cancel'])

    def test_normal_component_of_source_is_required_in_pressure(self):
        self.assertTrue(sourced_review()['checks']['normal_component_of_source_cannot_be_dropped'])

    def test_higher_harmonic_is_more_damped_on_the_same_domain(self):
        self.assertTrue(sourced_review()['checks']['higher_harmonic_damping_does_not_shrink_domain'])

    def test_liftwise_source_cannot_be_replaced_by_a_fiber_average(self):
        row=lifted_source_review()
        self.assertEqual(row['copies'],14)
        self.assertTrue(all(row['checks'].values()))
        self.assertGreater(row['midpoint_inverse_response_gap'],.1)

    def test_covering_counterexample_uses_the_actual_rectangle_not_a_shifted_interval(self):
        self.assertTrue(lifted_source_review()['checks']['every_path_starts_at_the_actual_rectangle_bottom'])

    def test_temporal_support_requires_an_actual_cutoff(self):
        row=failed_shortcuts_review()
        self.assertTrue(row['checks']['positive_finite_endpoint_is_not_zero_extension'])
        self.assertTrue(all(psi(u)==0 for u in (0,.05,.95,1)))
        self.assertAlmostEqual(psi(.5),1)

    def test_zero_harmonic_is_not_passed_to_the_wave_inverse(self):
        with mp.workdps(100):
            f=Fixture(51600,RADIUS,1)
            with self.assertRaises(ValueError):
                f.solve(m=0,source=True)


class JetsAndTailTests(unittest.TestCase):
    def test_all_fixed_homogeneous_derivative_degrees_are_finite(self):
        self.assertEqual(derivative_degrees(8),[2*j for j in range(9)])

    def test_source_Duhamel_cost_and_input_derivative_growth_are_retained(self):
        self.assertEqual(derivative_degrees(6,[0]*7),[2*j+1 for j in range(7)])
        self.assertEqual(derivative_degrees(6,[3*j for j in range(7)]),[3*j+1 for j in range(7)])

    def test_radius_dependent_tails_beat_every_fixed_physical_loss(self):
        self.assertTrue(all(derivative_and_tail_review()['checks'].values()))

    def test_zero_Gaussian_decay_cannot_prove_flatness(self):
        with self.assertRaises(ValueError):
            flat_tail_threshold(RADIUS,0,1,2,3,4)

    def test_old_radius_and_radius_uniform_inverse_shortcuts_fail(self):
        row=failed_shortcuts_review()
        self.assertTrue(all(row['checks'].values()))
        self.assertEqual(row['ideal_inverse_growth_when_radius_divided_by_four'],8)


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=generate()

    def test_all_required_gates_and_frozen_record(self):
        validate(self.record)
        frozen=json.loads((HERE/'evidence.json').read_text())
        validate(frozen)
        compare_records(self.record,frozen)

    def test_every_inherited_source_byte_is_preserved(self):
        hashes=inherited_provenance()
        self.assertGreater(len(hashes),404)
        self.assertTrue(all(hashes.values()))

    def test_scope_cannot_be_promoted_to_unconditional_blowup(self):
        row=copy.deepcopy(self.record)
        row['flags']['blowup_verified']=True
        with self.assertRaises(ValueError):
            validate(row)
        self.assertFalse(FLAGS['unconditional_complete_residual_improvement_verified'])

    def test_altered_radius_or_missing_source_support_condition_is_rejected(self):
        row=copy.deepcopy(self.record)
        row['protocol']['global_radius']='.04'
        with self.assertRaises(ValueError):
            validate(row)
        row=copy.deepcopy(self.record)
        del row['proof']['temporal_support']
        with self.assertRaises(ValueError):
            validate(row)

    def test_tiny_covariance_is_compared_relatively_not_erased_by_absolute_tolerance(self):
        row=copy.deepcopy(self.record)
        row['reviews']['full_Haar_primary_covariance']['rows'][0]['h']*=2
        with self.assertRaises(ValueError):
            compare_records(self.record,row)

    def test_existing_outputs_are_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)/'evidence.json'
            target.write_text('untouched')
            with self.assertRaises(SystemExit):
                main(['--output',str(target)])
            self.assertEqual(target.read_text(),'untouched')


if __name__=='__main__':
    unittest.main()
