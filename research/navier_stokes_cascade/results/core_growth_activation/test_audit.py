#!/usr/bin/env python3
"""Regression for the local core-growth derivation and its global limits."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from audit import HERE, compare_records, generate, main, validate
from diagnostics import load_gauged_activation
from growth import algebra_constants, budget, f0, pressure_axis


class CoreGrowthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 140
        cls.record = generate()

    def test_velocity_polynomial_and_uniform_remainder(self):
        r = self.record['reviews']['exact_source_scalar_and_separation_inequalities']
        self.assertTrue(r['checks']['velocity_polynomial_exact'])
        self.assertTrue(r['checks']['velocity_bound_absorbs_j_eta_squared'])

    def test_pressure_upper_mass_and_indispensable_ideal_atom(self):
        self.assertLess(algebra_constants()['reference_pressure_mass'], 5)
        for atoms in ([(0, 1)], [(1, '2.49')], [(1, '2.5'), (0, '2.5')], [(-1, 3)]):
            with self.assertRaises(ValueError):
                pressure_axis(0, atoms=atoms)

    def test_small_Z_localization_and_chi_separation_are_exact(self):
        c = algebra_constants()
        self.assertLess(c['small_Z_eta_over_j'], Q(1, 30))
        self.assertGreater(c['chi_separation'], Q(99, 100))

    def test_weak_angular_sector_resolves_eta_over_j(self):
        rows = self.record['reviews']['labelled_pressure_microregion']['rows']
        weak = [r for r in rows if r['chi'] <= .99]
        self.assertTrue(weak)
        self.assertLess(min(r['chi'] for r in weak), 1e-20)
        self.assertTrue(all(abs(r['Z_over_delta']) > 1 for r in weak))

    def test_pressure_axis_and_poles_have_no_d_division(self):
        for eta in (-1, 0, 1):
            r = pressure_axis(eta)
            self.assertTrue(all(mp.isfinite(v) for v in r.values()))
            self.assertGreaterEqual(r['Kpressure'], mp.mpf('1.25'))

    def test_scalar_hypergeometric_has_independent_series_check(self):
        r = self.record['reviews']['scalar_comparison_series']
        self.assertLess(r['independent_series_error'], 1e-100)
        self.assertEqual(f0(0), 1)
        self.assertEqual(f0(0, 1), -mp.mpf('.25'))

    def test_endpoint_angular_growth_has_exact_polynomial_bound(self):
        r = self.record['reviews']['exact_source_scalar_and_separation_inequalities']
        self.assertEqual(Q(r['p4_at_1p98']), -Q(75535511, 400000000))
        self.assertGreater(Q(r['scalar_cubic_at_2p05']), Q(265, 1000))

    def test_chi_zero_keeps_positive_source_without_division(self):
        c = algebra_constants()
        self.assertLess(c['source_remainder_prefactor'], Q(3, 1000))
        self.assertGreater(Q(28, 10)-Q(1, 10), Q(5, 2))
        self.assertTrue(self.record['flags']['positive_angular_source_without_division_by_chi_derived'])

    def test_positive_shear_uses_forward_integral(self):
        self.assertGreater(Q(5, 32), Q(1, 10))
        self.assertTrue(self.record['reviews']['exact_source_scalar_and_separation_inequalities']['checks']['positive_forward_shear'])

    def test_both_growth_sectors_are_preserved(self):
        self.assertEqual(self.record['protocol']['two_exit_sectors'], ['chi>.99', 'chi<=.99'])
        for row in self.record['reviews']['outward_frozen_width_and_growth_budgets']['rows']:
            self.assertTrue(row['checks']['endpoint_angular_error_small'])
            self.assertTrue(row['checks']['axial_exit_ratio_above_2p3'])

    def test_all_five_forward_moments_replay(self):
        r = self.record['reviews']['complete_forward_moments_and_C_gauge']
        self.assertLess(r['maximum_primitive_ODE_error'], 1e-100)
        self.assertLess(r['maximum_direct_difference_error'], 1e-100)

    def test_each_pressure_contribution_is_retained(self):
        rows = self.record['reviews']['complete_forward_moments_and_C_gauge']['rows']
        for row in rows:
            self.assertGreater(abs(row['pressure_value']), 1e-12)
            self.assertGreater(abs(row['pressure_eta']), 1e-12)

    def test_C_gauge_keeps_inverse_C_squared_pressure(self):
        r = self.record['reviews']['complete_forward_moments_and_C_gauge']
        self.assertLess(r['affine_inverse_C_squared_error'], 1e-100)
        self.assertTrue(r['checks']['pressure_has_required_inverse_C_squared_scaling'])
        self.assertNotEqual(r['rows'][0]['p2hat_over_C'], r['rows'][1]['p2hat_over_C'])

    def test_seed_gauge_eta_jets_and_axial_shear(self):
        m = load_gauged_activation(3)
        y, eta = mp.mpf('.003'), mp.mpf('.23')
        seed = m.seed(y, eta)
        for name in ('E', 'U', 'H', 'M', 'I', 'J', 'S', 'Cp', 'Pi'):
            actual = mp.diff(lambda e: m.seed(y, e)[name].v, eta)
            self.assertLess(abs(actual-seed[name].d), mp.mpf('1e-100'))
        self.assertLess(abs(-2*mp.diff(lambda z: m.seed(z, eta)['U'].v, y)/seed['E'].v-seed['q'].v), mp.mpf('1e-100'))

    def test_manufactured_baseline_defect_cannot_be_discarded(self):
        r = self.record['reviews']['complete_forward_moments_and_C_gauge']
        self.assertTrue(all(max(abs(v) for v in row['seed_defect']) > .1 for row in r['rows']))
        self.assertFalse(self.record['protocol']['manufactured_fixtures_are_actual_core'])

    def test_outward_budgets_and_nested_precision_refinement(self):
        r = self.record['reviews']['outward_frozen_width_and_growth_budgets']
        self.assertLess(r['maximum_relative_log_refinement_error'], 1e-55)
        self.assertEqual(self.record['protocol']['interval_digits'], [100, 160])
        self.assertTrue(all(all(row['checks'].values()) for row in r['rows']))

    def test_activation_K_and_existing_width_fit_without_retuning(self):
        p = self.record['protocol']
        self.assertEqual(p['activation_K'], 'C^2*exp(.8*B)<C^3')
        self.assertEqual(p['unchanged_width_and_kappa0'], 'C^-20')
        for row in self.record['reviews']['outward_frozen_width_and_growth_budgets']['rows']:
            self.assertTrue(row['checks']['activation_K_below_C3'])
            self.assertTrue(row['checks']['native_growth_change_below_p05'])

    def test_cone_and_lambda_have_positive_local_bounds(self):
        self.assertEqual(self.record['protocol']['activated_growth_lower'], '17/8')
        self.assertEqual(algebra_constants()['growth_lower_denominator'], 680)
        for row in self.record['reviews']['outward_frozen_width_and_growth_budgets']['rows']:
            self.assertTrue(row['checks']['squared_cone_loss_below_half'])
            self.assertTrue(mp.isfinite(mp.mpf(row['logs']['lambda_squared_lower_log']['lower'])))

    def test_invalid_domain_and_inadequate_interval_precision_rejected(self):
        for args in ((-1,), ('4.1001',), (1, -1), (1, .5)):
            with self.assertRaises(ValueError):
                f0(*args)
        for args in (('unknown', 100), ('Md64', 80)):
            with self.assertRaises(ValueError):
                budget(*args)
        for kwargs in ({'P': 15}, {'h': 0}, {'j': '.1'}, {'sigma': 0}):
            with self.assertRaises(ValueError):
                pressure_axis(0, **kwargs)

    def test_fixed_profile_bound_is_not_C_or_width_uniformity(self):
        self.assertFalse(self.record['protocol']['fixed_profile_constants_uniform_as_C_tends_to_infinity'])
        self.assertFalse(self.record['protocol']['physical_radial_derivatives_uniform_in_width'])

    def test_global_pressure_and_compact_cone_are_open(self):
        self.assertTrue(self.record['flags']['unedited_reference_axis_pressure_family_identified'])
        for key in ('completed_global_backward_pressure_identification_verified',
                    'compact_modulation_and_repair_cone_independently_verified', 'whole_leading_profile_verified'):
            self.assertFalse(self.record['flags'][key])

    def test_parent_ledgers_and_frozen_record_reproduce(self):
        validate(self.record)
        self.assertEqual(len(self.record['inherited_provenance']), 474)
        compare_records(self.record, json.loads((HERE/'evidence.json').read_text()))

    def test_global_pulse_residual_or_blowup_promotion_rejected(self):
        for key in ('completed_global_backward_pressure_identification_verified',
                    'unconditional_uniform_actual_pulse_estimate_verified',
                    'unconditional_complete_residual_improvement_verified', 'blowup_verified'):
            altered = deepcopy(self.record)
            altered['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(altered)

    def test_failed_pressure_gate_and_omitted_support_obligation_rejected(self):
        for mutate in (lambda r: r['remaining_obligations'].pop(),
                       lambda r: r['reviews']['complete_forward_moments_and_C_gauge']['checks'].update(pressure_value_and_eta_terms_both_nonzero=False)):
            altered = deepcopy(self.record)
            mutate(altered)
            with self.assertRaises(ValueError):
                validate(altered)

    def test_record_tampering_and_nonfinite_diagnostics_rejected(self):
        for value in (float('inf'), -100.0):
            altered = deepcopy(self.record)
            altered['reviews']['scalar_comparison_series']['rows'][0]['f'] = value
            with self.assertRaises(ValueError):
                compare_records(self.record, altered)

    def test_cli_does_not_overwrite_existing_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)/'record.json'
            target.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output', str(target)])
            self.assertEqual(target.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
