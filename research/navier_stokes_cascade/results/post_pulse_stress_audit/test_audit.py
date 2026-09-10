"""Post-pulse bounds, partial pressure and forbidden endpoint promotion."""
from copy import deepcopy
import json
import unittest
import mpmath as mp
from bounds import post_bound, HERE
from diagnostics import fixture, Patch, terminal_control, terminal_q_over_h, independent_terminal, Rule, relative
from audit import validate, compare_record, GATE_NAMES, STATUS


class BoundsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 280
        cls.b = post_bound()

    def test_all_constants(self):
        self.assertTrue(all(v is True for v in self.b['checks'].values()))

    def test_uniform_stress_bound_is_positive(self):
        self.assertTrue(0<self.b['w_universal_upper']<mp.mpf('1e-440'))
        self.assertGreater(self.b['second_expression_universal_upper'],0)

    def test_Q_floor_is_finite(self):
        self.assertEqual(self.b['Q_floor'],'h/2000')
        self.assertGreater(self.b['terminal_Q_over_h_lower'],mp.mpf(1)/2000)

    def test_positive_shear_is_separate(self):
        self.assertGreater(self.b['log_a_minus_two_lower_interval'][0],self.b['log_h_interval'][0])

    def test_finite_radius(self):
        self.assertGreater(self.b['universal_log_ps1_lower'],mp.log(100))
        self.assertGreater(self.b['finite_normalized_cone_gap_lower'],mp.mpf('1.84'))

    def test_edits_not_deleted(self):
        for k in ('angular_relative_edit_upper','angular_eta_log_edit_upper','angular_slope_edit_upper'):
            self.assertGreater(self.b[k],0)

    def test_unsupported_endpoint_and_inputs(self):
        for kw in ({'endpoint':'3'},{'endpoint':'infinity'},{'digits':20},{'radius':99},{'md':3}):
            with self.assertRaises(ValueError): post_bound(**kw)

    def test_scope_stays_local(self):
        self.assertIn('No heat collar',self.b['scope'])
        self.assertIn('global-profile-unverified',STATUS)


class PatchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 280
        cls.s = fixture(32)
        cls.patch = Patch(cls.s,'.5')
        cls.row = cls.patch.record([.02])

    def test_actual_moment_rows(self):
        self.assertLess(max(self.row['moment_relative_errors']),mp.mpf('1e-100'))

    def test_eta_derivative_system(self):
        self.assertLess(self.row['derivative_row_relative_error'],mp.mpf('1e-100'))
        self.assertTrue(any(v!=0 for v in self.row['coefficient_eta']))

    def test_eta_derivative_against_separate_coefficients(self):
        delta = mp.mpf('1e-20')
        plus = self.s.angular_correction(self.patch.eta+delta)['coefficients']
        minus = self.s.angular_correction(self.patch.eta-delta)['coefficients']
        for p,m,derivative in zip(plus,minus,self.row['coefficient_eta']):
            self.assertLess(relative((p-m)/(2*delta),derivative),mp.mpf('1e-36'))

    def test_partial_pressure_is_nonzero(self):
        self.assertGreater(abs(self.row['mid_patch_pressure_over_cscale']),mp.mpf('.001'))
        self.assertGreater(abs(self.row['gap_pressure_over_cscale']),mp.mpf('.001'))

    def test_pressure_restored_before_and_after(self):
        self.assertLess(abs(self.row['complete_pressure_relative_residual']),mp.mpf('1e-100'))
        self.assertEqual(self.row['post_patch_pressure'],0)

    def test_independent_pressure_ODE(self):
        self.assertLess(max(x['normalized_gap'] for x in self.row['independent_runs'][0]['rows']),mp.mpf('1e-9'))

    def test_quadratic_budget_positive(self):
        self.assertTrue(0<self.row['quadratic_normalized_positive_upper']<mp.mpf('1e-9'))

    def test_cached_discrepancy_matches_original(self):
        from diagnostics import Schedule
        eta = mp.mpf('.5')
        original = Schedule.interpolation_moment_discrepancy(self.s,eta)
        self.assertLess(relative(original,self.s.discrepancy_jet(eta)[0]),mp.mpf('1e-200'))
        with self.assertRaises(RuntimeError): self.s.pressure_parts(eta)


class TerminalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.t = terminal_control(320,48)

    def test_positive_until_exact_endpoint(self):
        self.assertTrue(all(x['Q_over_h']>0 for x in self.t['terminal_rows'][:-1]))
        self.assertEqual(self.t['terminal_rows'][-1]['Q_over_h'],0)

    def test_half_unit_floor(self):
        self.assertGreater(self.t['terminal_rows'][1]['Q_over_h'],mp.mpf(1)/2000)

    def test_backward_source_ODE(self):
        with mp.workdps(280):
            rows = independent_terminal(self.t['h'],Rule(32,8),.02)
            self.assertLess(max(x['normalized_gap'] for x in rows),mp.mpf('1e-9'))

    def test_default_terminal_co_uses_active_precision(self):
        with mp.workdps(100):
            rule = Rule(16,4)
            self.assertEqual(terminal_q_over_h(0,0,rule),terminal_q_over_h(0,0,rule,mp.mpf('.001')))

    def test_finite_h_difference_survives(self):
        self.assertLess(self.t['stable_finite_h_difference_over_h'],0)
        self.assertLess(self.t['direct_relative_gap'],mp.mpf('1e-40'))
        self.assertEqual(self.t['binary64_Qp_over_h_difference'],0)

    def test_tail_energy_term_cannot_be_dropped(self):
        self.assertLess(abs(self.t['pure_power_N_residual']),mp.mpf('1e-100'))
        self.assertGreater(self.t['pure_power_energy_term'],mp.mpf('.9'))
        self.assertGreater(abs(self.t['omitted_h_energy_N_over_E2']),mp.mpf('.9'))

    def test_infinite_reference_tail_fails_cone(self):
        self.assertEqual(self.t['tail_Pc'],0)
        self.assertFalse(self.t['tail_cone_pass'])


class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads((HERE/'evidence.json').read_text())

    def test_record_valid(self):
        validate(self.record)
        self.assertEqual(set(self.record['gates']),GATE_NAMES)

    def test_missing_gate(self):
        bad = deepcopy(self.record)
        bad['gates'].pop('infinite_tail_cone_failure')
        with self.assertRaises(ValueError): validate(bad)

    def test_endpoint_promotion(self):
        bad = deepcopy(self.record)
        bad['certificate']['terminal_endpoint'] = '3'
        with self.assertRaises(ValueError): validate(bad)

    def test_tiny_bound_zeroed(self):
        bad = deepcopy(self.record)
        bad['certificate']['w_universal_upper'] = '0'
        with self.assertRaises(ValueError): validate(bad)

    def test_partial_pressure_erased(self):
        bad = deepcopy(self.record)
        bad['patches'][0]['mid_patch_pressure_over_cscale'] = '0'
        with self.assertRaises(ValueError): validate(bad)

    def test_tail_failure_hidden(self):
        bad = deepcopy(self.record)
        bad['terminal_controls'][0]['tail_cone_pass'] = True
        with self.assertRaises(ValueError): validate(bad)

    def test_provenance_changed(self):
        bad = deepcopy(self.record)
        bad['provenance']['post_pulse_stress_audit/bounds.py'] = '0'*64
        with self.assertRaises(ValueError): validate(bad)

    def test_relative_comparison_of_tiny_terms(self):
        with mp.workdps(120):
            compare_record('1e-1000','1.00000000000000000001e-1000','/physical')
            with self.assertRaises(ValueError): compare_record('1e-1000','0','/physical')


if __name__=='__main__':
    unittest.main()
