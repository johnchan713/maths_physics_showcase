"""Meaningful equation, scale, failure-control and evidence checks."""
import importlib.util
import json
import math
import unittest
import mpmath as mp
from bounds import closure_bound, shape_enclosure, HERE
from pulse import Prefix, Continuation, MomentSolver, stable_two_row, family_increment, relative, Schedule
from audit import validate, compare_record, GATE_NAMES


class BoundsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=110
        cls.bound=closure_bound(512,60)

    def test_all_continuum_constants(self):
        self.assertTrue(all(self.bound['checks'].values()))

    def test_narrow_root_bracket(self):
        self.assertLess(self.bound['left_residual_upper'],0)
        self.assertGreater(self.bound['right_residual_lower'],0)
        self.assertGreater(self.bound['S_amplitude_derivative_lower'],mp.mpf('.35'))

    def test_nested_positive_shape_bounds(self):
        fine=shape_enclosure(1024,60)
        coarse=self.bound['shape']
        self.assertLess(coarse['K_lower'],fine['K_lower'])
        self.assertLess(fine['K_upper'],coarse['K_upper'])
        self.assertGreater(fine['prefix_lower'],0)
        self.assertGreater(fine['cutoff_lower'],0)

    def test_finite_errors_are_not_zero(self):
        self.assertGreater(self.bound['remainder_universal_upper'],0)
        self.assertLess(self.bound['remainder_universal_upper'],mp.mpf('1e-100'))
        self.assertGreater(self.bound['angular_contraction_upper'],0)

    def test_logarithmic_parameter_order(self):
        b=self.bound
        self.assertLess(b['log_h_interval'][1],b['log_lambda_interval'][0])
        self.assertLess(b['log_lambda_interval'][1],-256)

    def test_invalid_interval_budget(self):
        with self.assertRaises(ValueError): shape_enclosure(8,60)
        with self.assertRaises(ValueError): closure_bound(512,20)

    def test_no_global_profile_label(self):
        self.assertEqual(self.bound['status'],'reference-moments-closed-pulse-cone-unverified')


class ArithmeticTests(unittest.TestCase):
    def test_binary64_loses_distinct_rows(self):
        lam=1e-20
        self.assertEqual(math.exp(1-4*lam)-math.exp(1-2*lam),0)
        with mp.workdps(100):
            l=mp.mpf('1e-20')
            factor=mp.exp(1-2*l)*mp.expm1(-2*l)
            self.assertLess(factor,0)
            self.assertLess(abs(factor/l+2*mp.e),mp.mpf('1e-18'))

    def test_stable_solver_with_nearly_equal_slopes(self):
        with mp.workdps(160):
            lam=mp.mpf('1e-70')
            expected=[mp.mpf(2),mp.mpf('-.75')]
            masses=[mp.mpf('.3'),mp.mpf('.31')]
            slopes=[mp.mpf('.5')-lam,mp.mpf('.5')-2*lam]
            rhs=[mass*(expected[0]+mp.exp(2*s)*expected[1]) for mass,s in zip(masses,slopes)]
            actual=stable_two_row(lam,masses,rhs)
            self.assertLess(max(relative(a,b) for a,b in zip(actual,expected)),mp.mpf('1e-80'))

    def test_zero_lambda_rejected(self):
        with self.assertRaises(ValueError): stable_two_row(mp.mpf(0),[1,1],[1,2])

    def test_family_increment_precision_guard(self):
        with mp.workdps(100):
            with self.assertRaises(ValueError): family_increment()


class PulseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=100
        cls.prefix=Prefix('.0001')
        cls.solver=MomentSolver(cls.prefix)
        cls.post=Continuation(cls.prefix)
        cls.row=cls.solver.root('.5',cls.post)

    def test_both_moments_cancel(self):
        self.assertLess(max(self.row['linear_relative_errors']),mp.mpf('1e-80'))

    def test_full_S_quadratic(self):
        self.assertLess(self.row['quadratic_relative_residual'],mp.mpf('1e-80'))
        self.assertTrue(mp.mpf('.9')<self.row['amplitude']<mp.mpf('1.2'))

    def test_independent_source_integration(self):
        result=self.solver.source_check('.5',self.row['amplitude'])
        self.assertLess(result['maximum_error'],1e-9)

    def test_one_correction_is_insufficient(self):
        self.assertGreater(self.row['only_M_correction_J_relative_error'],mp.mpf('.8'))

    def test_leading_root_has_detectable_error(self):
        self.assertLess(self.row['leading_root_S_error'],mp.mpf('-.01'))
        self.assertGreater(self.row['amplitude']-self.row['leading_amplitude'],mp.mpf('.01'))

    def test_tiny_affine_and_energy_terms_are_retained(self):
        parts=self.row['S_parts']
        self.assertGreater(parts['bump_quadratic_increment'],0)
        self.assertNotEqual(parts['bump_linear_coefficient'],0)
        self.assertGreater(parts['bump_constant_increment'],0)
        self.assertTrue(all(v!=0 for v in self.row['affine_b_scaled']))

    def test_entire_future_energy_included(self):
        parts=self.row['post_mass']['parts']
        self.assertEqual(len(parts),8)
        self.assertEqual(parts[-1]['stage'],'infinite-power-tail')
        self.assertTrue(all(p['mass']>0 for p in parts))
        self.assertTrue(any(v!=0 for v in self.row['post_mass']['angular_energy_increments']))

    def test_cached_continuation_matches_original_calculation(self):
        expected=mp.mpf('722.5787625604125647781794062499000933879483069478468718599742875303845595240745887313356939187798648')
        self.assertLess(relative(self.row['post_mass']['total'],expected),mp.mpf('1e-70'))
        s=self.post.schedule
        for eta in (mp.mpf('.5'),mp.mpf(1)):
            original=Schedule.interpolation_moment_discrepancy(s,eta)
            cached=s.interpolation_moment_discrepancy(eta)
            self.assertLess(relative(original,cached),mp.mpf('1e-70'))
        with self.assertRaises(RuntimeError): s.pressure(0)

    def test_flat_prefix_omission_has_positive_bound(self):
        for entry in self.solver.main:
            self.assertTrue(0<entry['prefix_relative_error_upper']<mp.mpf('1e-100'))
            self.assertTrue(0<entry['cutoff_width']<entry['cutoff_saddle']<1)

    def test_full_amplitude_not_assumed_even(self):
        # The fixed swirl is even, but the affine axial prefix is odd.
        negative=self.solver.affine('-.5')
        positive=self.solver.affine('.5')
        self.assertEqual(negative['a'],positive['a'])
        for a,b in zip(negative['b'],positive['b']): self.assertEqual(a,-b)


class PrefixTests(unittest.TestCase):
    def test_against_previous_finite_parameter_reference(self):
        spec=importlib.util.spec_from_file_location('pulse_prefix_crosscheck',HERE.parent/'intermediate_decay_audit'/'reference.py')
        old=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(old)
        with mp.workdps(280):
            # This checks the moment identities at matched quadrature, not
            # convergence. The audit separately rejects this four-panel
            # resolution and requires sixteen panels for reported fixtures.
            new=Prefix(order=48,incoming_panels=4)
            previous=old.Reference(order=48)
            e=mp.mpf('.5')
            f=previous.power_factors(previous.Tw)
            E2=previous.state('power',previous.Tw,e)['E_squared']
            expected=dict(m0=e*f['K']/mp.sqrt(E2),j0=e*f['rK']/mp.sqrt(E2),
                          s0=e*e*f['K2']/E2-f['F']/2)
            actual=new.moments(e)
            self.assertLess(max(relative(actual[k],expected[k]) for k in expected),mp.mpf('1e-100'))
            self.assertLess(relative(mp.exp(2*new.log_eb)/(1+e*e)**2,E2),mp.mpf('1e-100'))


class RecordTests(unittest.TestCase):
    def setUp(self):
        self.record=json.loads((HERE/'evidence.json').read_text())

    def test_record_valid(self):
        validate(self.record)
        self.assertEqual(set(self.record['gates']),GATE_NAMES)

    def test_scope_cannot_be_promoted(self):
        self.record['status']='blowup-proved'
        with self.assertRaises(ValueError): validate(self.record)

    def test_gate_cannot_be_removed(self):
        self.record['gates'].pop('actual_fixture_M_J_rows')
        with self.assertRaises(ValueError): validate(self.record)

    def test_protocol_cannot_be_relaxed(self):
        self.record['protocol']['claimed_remainder_bound']='1'
        with self.assertRaises(ValueError): validate(self.record)

    def test_previous_failure_cannot_be_erased(self):
        self.record['old_lambda_failure']['lambda_w_squared']='0'
        with self.assertRaises(ValueError): validate(self.record)

    def test_positive_omission_bound_cannot_be_zeroed(self):
        self.record['family_diagnostic']['rows'][0]['omitted_bump_amplitude_error_upper']='0'
        with self.assertRaises(ValueError): validate(self.record)

    def test_incoming_resolution_failure_cannot_be_erased(self):
        self.record['controls']['incoming_resolution_failures'][0]['Ka_relative_gap']='0'
        with self.assertRaises(ValueError): validate(self.record)

    def test_nonfinite_evidence_rejected(self):
        self.record['fixtures'][0]['rows'][0]['amplitude']='NaN'
        with self.assertRaises(ValueError): validate(self.record)

    def test_comparison_retains_tiny_quantities(self):
        with self.assertRaises(ValueError):
            compare_record({'amplitude_increment':'1e-1000'},{'amplitude_increment':'0'})
        with self.assertRaises(ValueError):
            compare_record({'positive_bound':'1e-5000'},{'positive_bound':'2e-5000'})

    def test_log_scales_compare_in_original_relative_units(self):
        with self.assertRaises(ValueError):
            compare_record({'log_scale':'-1000000'},{'log_scale':'-1000000.0000000000001'})


if __name__=='__main__':
    unittest.main()
