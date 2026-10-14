#!/usr/bin/env python3
"""Band/phase/frame identities and safeguards around the conditional claim."""
import copy
from fractions import Fraction
import json
import unittest
import mpmath as mp
import numpy as np
from bounds import (barrier_certificate, cone_certificate, derivative_exponents,
                    endpoint_flat_threshold, eventual_band_threshold, exponent_ledger,
                    phase_majorant)
from charts import H, U, band, label, nearest_nonzero_integer, rational
from diagnostics import (band_frame_review, cone_and_edge_review, derivative_cutoff_review,
                         fixture_first_band, linear_remainder_review, partition_scaling_review,
                         phase_identity_review, pulse_covariance_review, support_review,
                         threshold_and_barrier_review)
from audit import FLAGS, HERE, compare_records, expected_protocol, main, proof_ledger, provenance, validate


class BandPhaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 120
        cls.phase, cls.frames = phase_identity_review(), band_frame_review()

    def test_implicit_normalization_through_third_order(self):
        self.assertTrue(self.phase['checks']['implicit_profile_to_chart_derivatives'])

    def test_actual_phase_gradient_and_transport_formula(self):
        self.assertTrue(self.phase['checks']['actual_phase_formula_chain_rule'])

    def test_fast_time_and_axial_phase_terms_cannot_be_dropped(self):
        self.assertTrue(self.phase['checks']['fast_time_term_required'])
        self.assertTrue(self.phase['checks']['axial_phase_derivative_required'])

    def test_full_moving_frame_and_left_inverse(self):
        self.assertTrue(self.frames['checks']['full_moving_frame_identity'])
        self.assertTrue(self.frames['checks']['analytic_moving_basis_derivative'])

    def test_phase_frame_scales_on_whole_pulse_fixture(self):
        self.assertTrue(self.frames['checks']['phase_and_frame_fixture_have_inverse_S_error'])
        self.assertTrue(self.frames['checks']['cover_and_carrier_scales'])

    def test_nonzero_rounding_bound_and_zero_guard(self):
        for x in map(mp.mpf, ('-2.1', '-.5', '-.1', '0', '.1', '.5', '2.1')):
            n = nearest_nonzero_integer(x)
            self.assertNotEqual(n, 0)
            self.assertLessEqual(abs(n-x), 1)
        self.assertTrue(self.frames['checks']['universal_half_integer_rounding_bound_rejected'])

    def test_moving_basis_term_is_retained(self):
        self.assertTrue(self.frames['checks']['moving_basis_derivative_required'])

    def test_high_precision_frame_values_reproduce(self):
        with mp.workdps(160):
            refined = band_frame_review()
        for a, b in zip(self.frames['rows'], refined['rows']):
            for name in ('scaled_normal_error', 'scaled_frame_entry_error', 'scaled_damping_error'):
                self.assertLess(abs(a[name]-b[name]), mp.mpf('1e-60'))


class InequalityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 120
        cls.thresholds = threshold_and_barrier_review()
        cls.cone = cone_and_edge_review()
        cls.derivatives = derivative_cutoff_review()

    def test_all_later_band_bound_has_analytic_threshold(self):
        self.assertTrue(self.thresholds['checks']['all_later_band_inequality'])
        self.assertEqual(fixture_first_band(), 51600)

    def test_low_bands_are_not_accepted(self):
        self.assertTrue(self.thresholds['checks']['low_bands_not_silently_accepted'])

    def test_invariant_interval_keeps_strict_radial_lower_bound(self):
        self.assertTrue(self.thresholds['checks']['riccati_barriers_point_inward'])
        self.assertTrue(self.thresholds['checks']['strictly_positive_radial_lower_bound'])

    def test_uncertified_barrier_and_cone_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            barrier_certificate('.6', 2, '.4', 1, '.3')
        with self.assertRaises(ValueError):
            cone_certificate('1.9', '1.5', '.2', 1)

    def test_nonfinite_or_invalid_threshold_inputs_are_rejected(self):
        for args in ((0, 1), ('.02', 1), ('.005', 0), ('.005', 'inf')):
            with self.assertRaises(ValueError):
                eventual_band_threshold(*args)
        with self.assertRaises(ValueError):
            endpoint_flat_threshold(1, 0, 1, 1, 1)

    def test_differentiated_duhamel_has_fixed_order_polynomial_losses(self):
        self.assertEqual(derivative_exponents(8), [2*j for j in range(9)])
        self.assertEqual(derivative_exponents(8, source_powers=[0]*9), [2*j+1 for j in range(9)])
        self.assertEqual(derivative_exponents(3, source_powers=[2, 5, 8, 11]), [3, 6, 9, 12])

    def test_all_fixed_endpoint_losses_are_absorbed(self):
        self.assertTrue(self.derivatives['checks']['all_fixed_endpoint_losses_absorbed'])

    def test_cone_neumann_bound_preserves_positive_weights(self):
        self.assertTrue(self.cone['checks']['continuum_cone_certificate'])
        self.assertTrue(self.cone['checks']['finite_cone_matrix_crosscheck'])
        self.assertTrue(self.cone['checks']['finite_neumann_error_crosscheck'])

    def test_target_outside_cone_is_not_promoted(self):
        self.assertTrue(self.cone['checks']['outside_cone_not_promoted_to_positive'])

    def test_square_root_edge_derivatives_retain_half_weight(self):
        self.assertTrue(self.cone['checks']['square_root_weight_derivatives_keep_half_weight'])
        self.assertTrue(self.cone['checks']['edge_weight_absorbs_every_fixed_inverse_power'])

    def test_exact_covariance_and_linear_gains_do_not_assert_full_residual_gain(self):
        ledger = exponent_ledger()
        self.assertEqual(ledger['real_primary_covariance_remainder'], Fraction(99999, 50000))
        self.assertEqual(ledger['conditional_linear_wave_gain'], Fraction(49997, 100000))
        self.assertGreater(ledger['relative_physical_stress_gain'], 0)
        self.assertLess(ledger['physical_covariance_remainder_power'], 0)
        self.assertFalse(ledger['full_residual_gain'])


class PulseAndOperatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 120
        cls.pulse = pulse_covariance_review()
        cls.operator = linear_remainder_review()
        cls.support = support_review()
        cls.partition = partition_scaling_review()

    def test_huge_band_stiff_pulse_has_positive_normalized_amplitude(self):
        for key in ('ratio_fits_supplied_entry_bound_interval', 'growing_amplitude_keeps_positive_lower_bound',
                    'gaussian_peak_resolved', 'actual_phase_fixture_covariance_directions'):
            self.assertTrue(self.pulse['checks'][key])

    def test_integrated_fixture_columns_have_positive_and_signed_inverses(self):
        self.assertTrue(self.pulse['checks']['integrated_primary_weights_positive'])
        self.assertTrue(self.pulse['checks']['integrated_signed_inverse_has_both_signs'])
        self.assertTrue(self.pulse['covariance_quadrature_is_truncated_diagnostic'])

    def test_tighter_stiff_solve_reproduces_covariance(self):
        refined = pulse_covariance_review(tighter=True)
        for a, b in zip(self.pulse['rows'], refined['rows']):
            np.testing.assert_allclose(a['covariance_column'], b['covariance_column'], rtol=2e-8, atol=2e-13)

    def test_finite_endpoint_pulse_is_small_but_not_zero(self):
        data = band(fixture_first_band())
        lab = label(data, 1)
        rho = U/2
        logp = lab['lambda0']*data['length']/U*(mp.asinh(rho)-mp.asinh(U)
            -((rho+rho**3/3)-(U+U**3/3))/(1+U*U)**mp.mpf('1.5'))
        self.assertLess(logp, -1000)
        self.assertGreater(mp.exp(logp), 0)
        self.assertTrue(proof_ledger()['global_linear_gain_requires_original_source_harmonic_smooth_temporal_extension'])

    def test_full_harmonic_linear_operator_matches_direct_differentiation(self):
        self.assertTrue(self.operator['checks']['full_principal_plus_linear_remainder_identity'])

    def test_mixed_angular_viscosity_and_pressure_omissions_fail(self):
        for key in ('angular_viscous_connection_retained', 'mixed_phase_amplitude_terms_retained',
                    'pressure_amplitude_gradient_retained'):
            self.assertTrue(self.operator['checks'][key])

    def test_common_path_and_cover_factor(self):
        self.assertTrue(self.support['checks']['pulse_path_fixes_slow_and_transverse_coordinates'])
        self.assertTrue(self.support['checks']['common_torus_cover_factor_required'])

    def test_distinct_preimages_must_not_be_averaged(self):
        self.assertTrue(self.support['checks']['allowed_sources_must_not_be_averaged_over_preimages'])

    def test_squared_partition_and_physical_normalization_are_retained(self):
        self.assertTrue(all(self.partition['checks'].values()))


class ScopeTests(unittest.TestCase):
    def test_protocol_and_all_inherited_bytes(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()), expected_protocol())
        inherited = provenance()
        self.assertEqual(len(inherited), 335)
        self.assertTrue(all(inherited.values()))

    def test_actual_global_field_and_full_nonlinear_cycle_remain_open(self):
        for key in ('actual_global_background_numerically_evaluated', 'actual_mean_correctors_constructed',
                    'full_nonlinear_interaction_and_residual_cycle_verified', 'smooth_force_verified', 'blowup_verified'):
            self.assertFalse(FLAGS[key])
        self.assertTrue(FLAGS['relies_on_inherited_continuum_background_leading_gaps_and_support_geometry'])

    def test_promotions_removed_gates_and_numeric_drift_are_rejected(self):
        record = json.loads((HERE/'evidence.json').read_text())
        validate(record)
        promoted = copy.deepcopy(record)
        promoted['flags']['blowup_verified'] = True
        with self.assertRaises(ValueError):
            validate(promoted)
        missing = copy.deepcopy(record)
        del missing['gates']['full_linear_operator']
        with self.assertRaises(ValueError):
            validate(missing)
        drift = copy.deepcopy(record)
        drift['reviews']['huge_band_pulse_fixture']['length'] *= 2
        with self.assertRaises(ValueError):
            compare_records(record, drift)

    def test_existing_evidence_is_not_overwritten(self):
        with self.assertRaises(SystemExit):
            main(['--output', str(HERE/'evidence.json')])


if __name__ == '__main__':
    unittest.main()
