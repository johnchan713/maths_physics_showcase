#!/usr/bin/env python3
"""Physical-coordinate, tensor, pulse and omitted-term crosschecks."""
import copy
import json
import math
import unittest
import numpy as np
from jets import Jet, arrays, bump, variables
from waves import advection, divergence, laplacian
from checks import (coordinate_review, covariance_review, residual_review,
                    mean_pressure_review)
from pulses import (gaussian_constants, log_envelope, pulse_review,
                    solve_pulse)
from audit import (FLAGS, HERE, REMAINING, compare_records, expected_protocol,
                   provenance, validate)


class TaylorAndCoordinateTests(unittest.TestCase):
    def test_cubic_mixed_derivatives_use_factorial_normalization(self):
        x, y, z, t = variables((.2, -.4, .7, .3))
        f = 2*x*x*y+3*x*y*z-4*z*t*t+x**3
        self.assertAlmostEqual(f.diff(0).diff(0).diff(1).value, 4.)
        self.assertAlmostEqual(f.diff(0).diff(1).diff(2).value, 3.)
        self.assertAlmostEqual(f.diff(2).diff(3).diff(3).value, -8.)
        self.assertAlmostEqual(f.diff(0).diff(0).diff(0).value, 6.)

    def test_analytic_composition_through_third_order(self):
        x, y, z, t = variables((.2, -.4, .7, .3))
        f = x.exp()*y.sin()+z.log()+t.atan()
        self.assertAlmostEqual(f.diff(0).diff(1).diff(1).value, -math.exp(.2)*math.sin(-.4))
        self.assertAlmostEqual(f.diff(2).diff(2).diff(2).value, 2/.7**3)
        self.assertAlmostEqual(f.diff(3).diff(3).diff(3).value, (6*.3**2-2)/(1+.3**2)**3)
        self.assertAlmostEqual(((1+x*x)**.5).diff(0).diff(0).value, (1+.2**2)**-1.5)

    def test_bump_is_zero_with_all_retained_boundary_and_exterior_jets(self):
        for x in (-1.5, -1., 1., 1.5):
            self.assertTrue(np.all(bump(Jet.variable(x, 0)).c == 0))

    def test_solid_rotation_vector_laplacian_and_advection(self):
        r, theta, z, t = variables((1.3, .2, -.1, .4))
        v = arrays([Jet(0), 2*r, Jet(0)])
        np.testing.assert_allclose(laplacian(v, r.value), np.zeros(3), atol=1e-14)
        np.testing.assert_allclose(advection(v, v, r.value), [-4*r.value, 0., 0.], atol=1e-14)
        self.assertAlmostEqual(divergence(v, r.value), 0.)


class PhysicalWaveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.coordinates = coordinate_review()
        cls.cov = covariance_review()
        cls.residual = residual_review()
        cls.mean = mean_pressure_review()

    def test_curl_matches_potential_and_cartesian_through_second_derivatives(self):
        self.assertTrue(self.coordinates['checks']['exact_curl_through_required_derivatives'])

    def test_complete_physical_momentum_matches_independent_coordinates(self):
        self.assertTrue(self.coordinates['checks']['full_physical_residual_independent_coordinates'])

    def test_transverse_amplitude_divergence_and_periodicity(self):
        for key in ('transverse_principal_amplitude', 'exact_divergence', 'global_angular_periodicity'):
            self.assertTrue(self.coordinates['checks'][key])

    def test_compact_physical_extension_away_from_axis(self):
        self.assertTrue(self.coordinates['checks']['smooth_compact_support_fixture'])

    def test_curl_connections_and_pressure_omissions_are_detected(self):
        for key in ('curl_remainder_cannot_be_discarded', 'radial_curl_connection_cannot_be_discarded',
                    'vector_laplacian_connections_cannot_be_discarded',
                    'advection_connections_cannot_be_discarded', 'pressure_cannot_be_discarded'):
            self.assertTrue(self.coordinates['checks'][key])

    def test_all_six_covariance_entries_match_angular_integrals(self):
        self.assertTrue(self.cov['checks']['full_six_entry_covariance'])

    def test_real_curl_covariance_error_is_quadratic_and_positive(self):
        self.assertTrue(self.cov['checks']['real_curl_covariance_error_is_quadratic'])
        self.assertTrue(self.cov['checks']['curl_covariance_error_positive_semidefinite'])

    def test_full_tensor_divergence_includes_cylindrical_connections(self):
        self.assertTrue(self.cov['checks']['full_tensor_mean_divergence'])

    def test_signed_inverse_retains_every_full_velocity_increment_term(self):
        self.assertTrue(self.cov['checks']['signed_inverse_and_all_covariance_increment_terms'])
        self.assertTrue(self.cov['checks']['signed_linear_inverse_does_not_remove_quadratic_remainder'])

    def test_distinct_phases_do_not_imply_zero_mean_cross_terms(self):
        self.assertTrue(self.cov['checks']['angular_mean_cross_terms_require_a_reason'])

    def test_full_residual_retains_wave_background_and_wave_wave_terms(self):
        self.assertTrue(self.residual['checks']['full_residual_includes_all_wave_background_and_wave_wave_terms'])
        self.assertTrue(self.residual['checks']['full_mean_residual_is_full_tensor_divergence'])

    def test_fixed_covariance_kinematic_field_has_uncancelled_viscosity(self):
        for key in ('kinematic_covariance_matching_does_not_reduce_residual',
                    'uncancelled_viscosity_has_k_squared_scale',
                    'finite_frequency_values_approach_derived_limit'):
            self.assertTrue(self.residual['checks'][key])
        cross = np.array([r['principal_cross_covariance'] for r in self.residual['rows']])
        np.testing.assert_allclose(cross, np.broadcast_to(cross[0], cross.shape), rtol=0., atol=1e-14)

    def test_radial_pressure_does_not_cancel_the_axial_mean_equation(self):
        self.assertTrue(self.mean['checks']['radial_pressure_leaves_axial_defect'])

    def test_radial_pressure_primitive_can_leave_inner_pressure_constant(self):
        self.assertTrue(self.mean['checks']['pressure_primitive_need_not_vanish_inside_wave_support'])

    def test_tangential_stress_matching_leaves_other_mean_fluxes(self):
        self.assertTrue(self.mean['checks']['tangential_flux_matching_leaves_full_mean_terms'])
        self.assertTrue(self.mean['checks']['theta_z_flux_cannot_be_discarded'])


class PulseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = pulse_review()

    def test_moving_normal_constraint_pressure_and_energy(self):
        for key in ('moving_normal_constraint', 'complex_pressure_principal_equation',
                    'shear_viscosity_energy_identity'):
            self.assertTrue(self.record['checks'][key])

    def test_closed_envelope_gaussian_and_midpoint(self):
        self.assertTrue(self.record['checks']['continuum_gaussian_bound_sample_crosscheck'])
        c, bigc = gaussian_constants()
        self.assertGreater(c, 0.)
        self.assertGreater(bigc, c)
        for length in (40, 100, 240):
            self.assertAlmostEqual(log_envelope(length/2, length), 0.)
            self.assertLess(log_envelope(0., length), 0.)
            self.assertLess(log_envelope(length, length), 0.)

    def test_integrated_columns_have_positive_primary_weights_and_reference_scale(self):
        for key in ('positive_covariance_inverse_in_reference_fixture',
                    'covariance_directions_approach_reference', 'weight_has_reference_sqrt_length_scale'):
            self.assertTrue(self.record['checks'][key])

    def test_integrated_columns_admit_signed_linear_change(self):
        self.assertTrue(self.record['checks']['signed_linear_inverse_retained'])
        self.assertLess(max(r['signed_inverse_error'] for r in self.record['rows']), 1e-11)

    def test_tighter_solve_reproduces_covariance_columns(self):
        refined = pulse_review(tighter=True)
        for old, new in zip(self.record['rows'], refined['rows']):
            np.testing.assert_allclose(old['covariance_columns'], new['covariance_columns'],
                                       rtol=2e-9, atol=2e-11)

    def test_missing_moving_normal_pressure_i_and_viscosity_fail(self):
        for key in ('moving_normal_omission_detected', 'pressure_i_omission_detected',
                    'viscous_damping_omission_detected'):
            self.assertTrue(self.record['checks'][key])

    def test_normalized_solver_does_not_start_with_tiny_physical_initial_data(self):
        solution, _ = solve_pulse(240, 1)
        self.assertGreater(np.linalg.norm(solution.y[:, 0]), 1.)
        self.assertLess(math.exp(log_envelope(0., 240)), 1e-10)
        self.assertFalse(self.record['uniform_actual_pulse_or_parameter_derivative_bound_verified'])


class ScopeAndProvenanceTests(unittest.TestCase):
    def test_actual_pulses_mean_full_residual_and_blowup_are_open(self):
        for key in ('actual_band_phase_and_pulse_uniform_bounds_verified', 'actual_mean_correctors_constructed',
                    'full_physical_residual_improvement_verified', 'smooth_force_verified', 'blowup_verified'):
            self.assertFalse(FLAGS[key])
        self.assertIn('actual_mean_corrections_integral_constraints_and_all_interactions', REMAINING)

    def test_protocol_is_frozen(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())

    def test_all_inherited_audit_bytes_are_preserved(self):
        checks = provenance()
        self.assertGreaterEqual(len(checks), 320)
        self.assertTrue(all(checks.values()))

    def test_evidence_rejects_promotion_and_numerical_drift(self):
        record = json.loads((HERE/'evidence.json').read_text())
        validate(record)
        promoted = copy.deepcopy(record)
        promoted['flags']['blowup_verified'] = True
        with self.assertRaises(ValueError):
            validate(promoted)
        drifted = copy.deepcopy(record)
        drifted['reviews']['full_residual']['rows'][0]['full_increment_rms'] *= 2
        with self.assertRaises(ValueError):
            compare_records(record, drifted)


if __name__ == '__main__':
    unittest.main()
