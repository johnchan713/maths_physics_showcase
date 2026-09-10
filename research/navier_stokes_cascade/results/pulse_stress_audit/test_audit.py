"""Continuum assumptions, independent equations and evidence failure controls."""
from copy import deepcopy
from fractions import Fraction
import json
import math
import unittest
import mpmath as mp
from bounds import pulse_bound, HERE
from diagnostics import (SyntheticPulse, shape_jet, step_jet, algebra_controls,
                         affine_control, finite_cone_controls, angular_equality_control)
from audit import validate, compare_record, GATE_NAMES, STATUS


class BoundsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 120
        cls.bound = pulse_bound()

    def test_every_outward_constant(self):
        self.assertTrue(all(v is True for v in self.bound['checks'].values()))

    def test_strict_ratio_margins(self):
        self.assertGreater(self.bound['first_A24_margin_lower'],mp.mpf('1.48'))
        self.assertGreater(self.bound['second_A24_margin_lower'],mp.mpf('.82'))

    def test_positive_finite_remainders(self):
        b = self.bound
        self.assertTrue(all(v>0 for v in b['w_remainder_parts'].values()))
        self.assertGreater(b['w_remainder_universal_upper'],0)
        self.assertLess(b['w_remainder_universal_upper'],mp.mpf('1e-48'))

    def test_finite_radius_is_checked(self):
        self.assertGreater(self.bound['universal_log_ps1_lower'],mp.log(20000000))
        self.assertGreater(self.bound['finite_normalized_quadratic_gap_lower'],mp.mpf('.55'))

    def test_logarithmic_finite_parameter_order(self):
        b = self.bound
        self.assertLess(b['log_h_interval'][1],b['log_lambda_interval'][0])
        self.assertGreater(b['log_v_minus_two_lower_interval'][0],b['log_lambda_interval'][0])

    def test_pulse_scope_does_not_repair_Md4(self):
        b = pulse_bound(md=4)
        self.assertEqual(b['status'],'corrected-reference-pulse-cone-bounded')
        self.assertIn('not certified',b['scope'])
        old = json.loads((HERE.parent/'axial_stress_audit'/'evidence.json').read_text())
        self.assertLess(mp.mpf(old['old_Md_failure']['Pc_over_ps1']),0)

    def test_small_angle_inequality(self):
        # eta=t*sqrt(lambda) reduces the comparison to exact rationals.
        # No floating <= test is allowed at its exact equality point.
        for scale in ('0','0.1','1','3','1000'):
            t = Fraction(scale)
            self.assertEqual(1+t*t-2*t,(t-1)**2)
            self.assertLessEqual(t/(1+t*t),Fraction(1,2))

    def test_invalid_inputs_rejected(self):
        for kwargs in ({'md':3},{'md':4.5},{'digits':30},{'radius':99}):
            with self.assertRaises(ValueError): pulse_bound(**kwargs)


class AlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = algebra_controls()
        cls.cone = finite_cone_controls()

    def test_full_source_identity(self):
        self.assertLess(self.data['maximum_source_identity_error'],mp.mpf('1e-80'))

    def test_finite_h_cannot_be_dropped(self):
        self.assertGreater(self.data['omitted_h_maximum_gap'],mp.mpf('1e-6'))

    def test_angular_derivative_cannot_be_dropped(self):
        self.assertGreater(self.data['omitted_eta_maximum_gap'],mp.mpf('1e-6'))

    def test_energy_cannot_be_dropped(self):
        self.assertGreater(self.data['omitted_energy_maximum_gap'],mp.mpf('1e-6'))

    def test_finite_cone_algebra(self):
        self.assertLess(self.cone['identity_maximum_error'],mp.mpf('1e-80'))

    def test_ratios_alone_are_insufficient(self):
        small = self.cone['small_radius']
        self.assertTrue(small['passes_ratio_test'])
        self.assertFalse(small['passes_finite_cone'])
        self.assertLess(small['Pc'],small['v'])


class DynamicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.profile = SyntheticPulse(.3)
        cls.coarse = cls.profile.integrate()
        cls.fine = cls.profile.integrate(1,.02)

    def test_independent_source_ODE(self):
        self.assertLess(max(x['normalized_gap'] for x in self.coarse+self.fine),1e-9)

    def test_step_refinement(self):
        for x,y in zip(self.coarse,self.fine):
            for key in ('Q_source','N_over_E_source'):
                self.assertLess(abs(x[key]-y[key])/(1+abs(y[key])),1e-9)

    def test_actual_signed_patches_exercised(self):
        self.assertGreater(max(x['omitted_bump_bs_gap'] for x in self.fine),1e-6)
        first,second = [self.profile.profile(y)[1] for y in self.profile.centers]
        self.assertGreater(first,0)
        self.assertLess(second,0)

    def test_smooth_shape_support(self):
        for x in (-1,0,11,12): self.assertEqual(shape_jet(x),(0.0,0.0))
        self.assertEqual(step_jet(0),(0.0,0.0,0.0))
        self.assertEqual(step_jet(1),(1.0,0.0,0.0))

    def test_cutoff_derivative_sign(self):
        _,derivative = shape_jet(10.5)
        self.assertLess(derivative,-70)
        for xi in (.001,.005,.01,.015,.02,.5,1,10,10.3,10.5,10.9):
            R,Rx = shape_jet(xi)
            self.assertGreaterEqual(R,0)
            self.assertLessEqual(Rx,1)

    def test_non_even_strength_is_not_silently_symmetrized(self):
        self.assertNotEqual(SyntheticPulse(.3).profile(25)[1],SyntheticPulse(-.3).profile(25)[1])


class ArithmeticTests(unittest.TestCase):
    def test_80_digits_can_lose_the_stress(self):
        low = affine_control(80)
        self.assertEqual(low['direct_w'],0)
        self.assertLess(low['stable_w'],-1)
        self.assertGreater(low['direct_relative_error'],mp.mpf('.9'))

    def test_180_and_260_digits_recover_the_cancellation(self):
        for digits in (180,260):
            for scale in ('0','.5','2','10'):
                row = affine_control(digits,scale)
                self.assertLess(row['direct_relative_error'],mp.mpf('1e-60'))

    def test_first_convolution_derivative_matters(self):
        self.assertGreater(affine_control(260)['omitted_derivative_absolute_error'],2)

    def test_strict_excess_must_be_separated(self):
        row = affine_control(80)
        self.assertGreater(row['a_minus_two'],0)
        self.assertEqual(row['rounded_a_minus_two'],0)
        self.assertEqual(row['binary64_rounded_a_minus_two'],0)

    def test_actual_Md64_binary64_underflow(self):
        loglam = pulse_bound()['log_lambda_interval'][1]
        self.assertTrue(mp.isfinite(loglam))
        self.assertEqual(math.exp(float(loglam)),0)

    def test_false_floating_equality_comparison_is_retained(self):
        row = angular_equality_control()
        self.assertTrue(row['exact_rational_checks'])
        self.assertFalse(row['floating_comparison_passed'])
        self.assertTrue(0<row['relative_rounding_excess']<mp.mpf('1e-110'))


class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads((HERE/'evidence.json').read_text())

    def test_frozen_record(self):
        validate(self.record)
        self.assertEqual(set(self.record['gates']),GATE_NAMES)

    def test_missing_gate_rejected(self):
        bad = deepcopy(self.record)
        bad['gates'].pop('finite_radius_cone')
        with self.assertRaises(ValueError): validate(bad)

    def test_scope_promotion_rejected(self):
        bad = deepcopy(self.record)
        bad['status'] = 'global-blowup-proved'
        with self.assertRaises(ValueError): validate(bad)
        self.assertIn('global-profile-unverified',STATUS)

    def test_code_hash_change_rejected(self):
        bad = deepcopy(self.record)
        bad['provenance']['pulse_stress_audit/bounds.py'] = '0'*64
        with self.assertRaises(ValueError): validate(bad)

    def test_changed_margin_with_true_gate_rejected(self):
        bad = deepcopy(self.record)
        bad['certificate']['first_A24_margin_lower'] = '-1'
        with self.assertRaises(ValueError): validate(bad)

    def test_erased_tiny_error_with_true_gate_rejected(self):
        bad = deepcopy(self.record)
        bad['certificate']['w_remainder_parts']['end_corrections'] = '0'
        with self.assertRaises(ValueError): validate(bad)

    def test_changed_ODE_measurement_with_true_gate_rejected(self):
        bad = deepcopy(self.record)
        bad['ode_runs'][0]['rows'][0]['normalized_gap'] = .1
        with self.assertRaises(ValueError): validate(bad)

    def test_inherited_failure_cannot_be_erased(self):
        bad = deepcopy(self.record)
        bad['inherited']['old_Md_failure']['Pc_over_ps1'] = '1'
        with self.assertRaises(ValueError): validate(bad)

    def test_relative_comparison_preserves_tiny_physical_terms(self):
        with mp.workdps(100):
            compare_record('1e-1000','1.00000000000000000001e-1000','/physical')
            with self.assertRaises(ValueError): compare_record('1e-1000','0','/physical')
            with self.assertRaises(ValueError): compare_record('1e-1000','1.01e-1000','/physical')


if __name__=='__main__':
    unittest.main()
