#!/usr/bin/env python3
"""Direct heat/viscosity signs, edge ratios and scope regression checks."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from audit import FLAGS, HERE, compare_records, generate, main, validate
from tail import HeatTail, axial_coefficients, heat_derivative, positive_coefficient, step, step_jets


class OuterEdgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.record = generate()
        cls.tail = HeatTail(order=24, panels=8)

    def test_exact_coefficient_cancellation(self):
        self.assertTrue(all(self.record['reviews']['exact_axial_coefficient_algebra']['checks'].values()))

    def test_interval_bounds_and_actual_schedule_amplitude(self):
        result = self.record['reviews']['whole_interval_heat_and_schedule_bounds']
        self.assertTrue(all(result['checks'].values()))
        self.assertLess(result['positive_integral_coefficient_upper'], 7)
        self.assertLess(result['log_e_star_over_h14_upper'], -469)

    def test_exact_heat_integral_and_heat_equation(self):
        self.assertTrue(all(self.record['reviews']['exact_heat_integral']['checks'].values()))
        self.assertEqual(heat_derivative('.005', 0), 1)
        self.assertEqual(heat_derivative('.005', 0, 1), -mp.mpf('.005')*mp.mpf('1.005'))

    def test_heat_derivative_bounds_from_gamma_moments(self):
        h, Z = mp.mpf('.005'), mp.mpf('.02')
        for j in range(1, 6):
            value = heat_derivative(h, Z, j)
            self.assertGreater((-1)**j*value, 0)
            self.assertLessEqual(abs(value), mp.rf(h, j)*mp.rf(1+h, j))

    def test_flat_step_derivatives_against_independent_differentiation(self):
        for x in map(mp.mpf, ('.15', '.4', '.5', '.7', '.95')):
            _, first, second = step_jets(x)
            self.assertLess(abs(first-mp.diff(step, x)), mp.mpf('1e-65'))
            self.assertLess(abs(second-mp.diff(step, x, 2)), mp.mpf('1e-65'))

    def test_signed_source_and_positive_stress_agree(self):
        result = self.record['reviews']['signed_and_positive_viscous_stress']
        self.assertTrue(all(result['checks'].values()))
        for row in result['rows']:
            self.assertGreater(row['T1positive'], 0)

    def test_nonunit_q_physical_source_and_cylindrical_sign(self):
        result = self.record['reviews']['physical_axial_source_and_radial_stress_sign']
        self.assertTrue(all(result['checks'].values()))
        self.assertEqual([r['q'] for r in result['rows']], [1., .4, .2, .7])

    def test_beta_cannot_be_dropped_at_equator(self):
        _, alpha, beta = axial_coefficients('.005', 0)
        self.assertEqual(alpha, 0)
        self.assertEqual(beta, -2)
        row = self.record['reviews']['physical_axial_source_and_radial_stress_sign']['rows'][1]
        self.assertGreater(row['omitted_beta_source_gap'], 1)

    def test_positive_integral_coefficient_on_both_parameter_boundaries(self):
        for eta in (-1, 0, 1):
            self.assertGreaterEqual(positive_coefficient('.005', eta, 0), 2)

    def test_cone_loss_decreases_for_nonnegative_corrections(self):
        for row in self.record['reviews']['signed_and_positive_viscous_stress']['rows']:
            losses = [r['directional_loss'] for r in row['updated']]
            self.assertTrue(all(a >= b for a, b in zip(losses, losses[1:])))

    def test_flat_stress_ratio_uses_nonzero_coefficients_after_float_underflow(self):
        result = self.record['reviews']['flat_edge_and_correction_ratio']
        self.assertTrue(all(result['checks'].values()))
        for row in result['rows']:
            if row['delta'] == '1e-7':
                self.assertEqual(row['binary64_flat_factor'], 0)
                self.assertGreater(row['b0'], 0)

    def test_endpoint_ratio_retains_the_large_Xb_factor(self):
        eta = mp.mpf('.5')
        edge = self.tail.endpoint(eta)
        L = 1-2*self.tail.h*eta*eta
        wrong = 2*eta*eta/L**2
        self.assertEqual(edge['correction_ratio']/wrong, self.tail.xb)
        self.assertGreater(edge['correction_ratio'], 10000)

    def test_equatorial_endpoint_ratio_zero_but_interior_correction_positive(self):
        self.assertEqual(self.tail.endpoint(0)['correction_ratio'], 0)
        self.assertGreater(self.tail.factored('.001', 0)['b1'], 0)

    def test_raw_zero_edge_ratio_is_rejected(self):
        with self.assertRaises(ValueError):
            self.tail.direct(3, 0)
        with self.assertRaises(ValueError):
            self.tail.factored(0, 0)
        self.assertEqual(self.tail.endpoint(1)['stress'], [0, 0])
        self.assertEqual(self.tail.endpoint(1)['direction'], [1, 0])

    def test_nonnegative_Borel_cutoff(self):
        for s in map(mp.mpf, ('0', '.25', '.5', '.6', '.9', '1', '2')):
            chi = 1-step(2*s-1)
            self.assertGreaterEqual(chi, 0)
            self.assertLessEqual(chi, 1)

    def test_weighted_smallness_and_negative_cutoff_shortcuts_fail(self):
        result = self.record['reviews']['failed_shortcuts']
        self.assertTrue(all(result['checks'].values()))
        self.assertLess(result['generic_signed_perturbation_over_leading'], 0)

    def test_refinement_preserves_tiny_stress_normalization(self):
        self.assertTrue(self.record['refinement']['check'])
        altered = deepcopy(self.record)
        altered['reviews']['signed_and_positive_viscous_stress']['rows'][0]['T1positive'] = 0.
        with self.assertRaises(ValueError):
            compare_records(self.record, altered)

    def test_invalid_kernel_and_geometry_inputs_rejected(self):
        for h, z, j in ((0, 0, 0), ('.01', 0, 0), ('.005', -1, 0), ('.005', 0, -1)):
            with self.assertRaises(ValueError):
                heat_derivative(h, z, j)
        with self.assertRaises(ValueError):
            HeatTail(xtail=99)
        with self.assertRaises(ValueError):
            self.tail.state(2, '1.01')

    def test_all_frozen_hashes_and_new_record_reproduce(self):
        validate(self.record)
        self.assertTrue(all(self.record['inherited_provenance'].values()))
        compare_records(self.record, json.loads((HERE/'evidence.json').read_text()))

    def test_global_identification_or_blowup_promotion_rejected(self):
        for key in ('global_forward_backward_moment_identification_independently_reverified',
                    'inner_core_collar_and_compact_cone_foundation_independently_verified',
                    'unconditional_uniform_pulse_estimate_verified', 'blowup_verified'):
            altered = deepcopy(self.record)
            altered['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(altered)

    def test_omitted_gate_or_remaining_obligation_rejected(self):
        altered = deepcopy(self.record)
        altered['gates'].pop('physical_axial_source_and_radial_stress_sign')
        with self.assertRaises(ValueError):
            validate(altered)
        altered = deepcopy(self.record)
        altered['remaining_obligations'].pop()
        with self.assertRaises(ValueError):
            validate(altered)

    def test_cli_never_overwrites_existing_record(self):
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp)/'record.json'
            target.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output', str(target)])
            self.assertEqual(target.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
