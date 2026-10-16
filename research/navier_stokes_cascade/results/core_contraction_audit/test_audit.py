#!/usr/bin/env python3
"""Core contraction, original equations, exact stress and scope regression."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from audit import HERE, compare_records, generate, main, validate
from core import Series, core_budget, finite_angular_inverse, formal_core, mixed_allocation, source_sums, weight


class CoreContractionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 140
        cls.record = generate()

    def test_exact_coefficient_shift_ratios(self):
        self.assertTrue(all(self.record['reviews']['coefficient_weights_and_derivative_allocation']['checks'].values()))

    def test_mixed_derivative_at_boundary_indices(self):
        for nu in (1, 2):
            self.assertEqual(mixed_allocation(0, 0, nu), 0)
            self.assertEqual(mixed_allocation(10**6, 0, nu), 0)
            self.assertLessEqual(mixed_allocation(0, 10**6, nu), 1)

    def test_weight_rho_dependence_is_preserved(self):
        self.assertEqual(weight(3, 4, Q(1, 7))/weight(3, 4), 7**4)

    def test_complex_poles_and_pressure_family_bounds(self):
        self.assertTrue(all(self.record['reviews']['pressure_family_complex_neighborhood']['checks'].values()))

    def test_complete_source_budgets_and_pressure_constant(self):
        r = source_sums()
        self.assertEqual(r['R1_linear'], Q(19237, 50))
        self.assertEqual(r['R1_quadratic_without_zeta'], Q(952, 25))
        self.assertEqual(r['R2_pressure'], 34)
        self.assertTrue(all(self.record['reviews']['complete_source_and_evaluation_budgets']['checks'].values()))

    def test_radial_inverses_solve_regular_equations(self):
        Series.ny, Series.nk = 8, 4
        F = Series({(0, 0): Q(3, 7), (3, 2): Q(5, 11)})
        Y = Series({(1, 0): 1})
        for nu in (1, 2):
            G = F.J(nu)
            self.assertTrue((Y*G.radial().radial()+nu*G.radial()-F).zero(6, 4))
            self.assertFalse(any(n == 0 for n, _ in G.terms))

    def test_angular_inverse_is_two_sided_without_small_chi(self):
        self.assertTrue(all(self.record['reviews']['order_one_angular_inverse']['checks'].values()))

    def test_angular_inverse_cannot_be_replaced_by_identity(self):
        Series.ny, Series.nk = 4, 2
        one, chi = Series.cast(1), Series.cast(3)
        inverse = finite_angular_inverse(one, chi)
        self.assertNotEqual(inverse.terms.get((1, 0)), 0)
        self.assertTrue((inverse+(chi*inverse).J(2)/2-one).zero(4, 2))

    def test_original_core_equations_cancel_with_eta_padding(self):
        r = self.record['reviews']['original_core_equations_and_forward_stress_jets']
        self.assertTrue(r['checks']['original_equations_and_forward_stress_cancel_in_interior_jets'])
        self.assertTrue(all(row['eta_padding'] >= row['radial_degree'] for row in r['rows']))

    def test_forward_primitives_cancel_both_stress_components(self):
        r = formal_core(3, 2)
        self.assertTrue(r['angular_stress'].zero(2, 2))
        self.assertTrue(r['axial_stress_times_E'].zero(2, 2))

    def test_pressure_and_azimuthal_datum_kept(self):
        r = self.record['reviews']['original_core_equations_and_forward_stress_jets']
        self.assertTrue(r['checks']['amplitude_obeys_azimuthal_datum'])
        self.assertTrue(all(row['identities']['pressure'] for row in r['rows']))

    def test_each_omitted_pressure_group_breaks_original_equations(self):
        r = self.record['reviews']['original_core_equations_and_forward_stress_jets']['omitted_pressure_controls']
        self.assertEqual(set(r), {'zero', 'drop_eta', 'drop_radial'})
        self.assertTrue(all(all(v.values()) for v in r.values()))

    def test_regular_axis_cannot_hide_nonzero_integration_constants(self):
        Series.ny, Series.nk = 4, 2
        with self.assertRaises(ValueError):
            Series({(0, 0): 1}).shift(-1)
        with self.assertRaises(ValueError):
            Series({(1, 0): 1}).shift(-2)

    def test_universal_contraction_ceiling_is_not_a_parameter_scan(self):
        r = self.record['reviews']['complete_source_and_evaluation_budgets']
        self.assertEqual(Q(r['universal_error_ceiling']), Q(1, 10**50))
        self.assertEqual(Q(8)-Q(1, 16), Q(127, 16))
        self.assertTrue(self.record['protocol']['finite_diagnostics_are_not_all_index_or_existence_proofs'])

    def test_independent_interval_budgets_match_selected_parameters(self):
        r = self.record['reviews']['independent_parameter_arithmetic']
        self.assertLess(r['maximum_relative_log_refinement_error'], 1e-55)
        self.assertTrue(all(r['checks'].values()))

    def test_core_regularities_do_not_claim_growth_or_full_residual(self):
        for name in ('actual_core_positive_shear_growth_and_activation_K_verified',
                     'whole_leading_profile_verified', 'unconditional_complete_residual_improvement_verified'):
            self.assertFalse(self.record['flags'][name])

    def test_real_axis_supremum_cannot_replace_analytic_norm(self):
        r = self.record['reviews']['failed_real_only_and_full_residual_shortcuts']
        self.assertEqual(r['real_supremum_one_analytic_norm_lower_bounds'], [4, 40, 400])

    def test_finite_formal_jets_are_not_actual_giant_schedule(self):
        r = self.record['reviews']['original_core_equations_and_forward_stress_jets']
        self.assertFalse(r['actual_giant_schedule_solution_evaluated'])
        self.assertFalse(self.record['protocol']['formal_diagnostics_are_actual_core'])

    def test_invalid_indices_scales_and_inverse_inputs_rejected(self):
        for n, k in ((-1, 0), (0, -1), (Q(1, 2), 0)):
            with self.assertRaises(ValueError):
                weight(n, k)
        for args in (('unknown', 80), ('P16', 10)):
            with self.assertRaises(ValueError):
                core_budget(*args)
        with self.assertRaises(ValueError):
            Series().inverse()
        with self.assertRaises(ValueError):
            formal_core(0, 2)

    def test_parent_ledgers_and_frozen_record_reproduce(self):
        validate(self.record)
        self.assertTrue(all(self.record['inherited_provenance'].values()))
        compare_records(self.record, json.loads((HERE/'evidence.json').read_text()))

    def test_pressure_identification_growth_and_blowup_promotion_rejected(self):
        for key in ('actual_global_pressure_family_identification_independently_verified',
                    'actual_core_positive_shear_growth_and_activation_K_verified',
                    'unconditional_uniform_actual_pulse_estimate_verified', 'blowup_verified'):
            altered = deepcopy(self.record)
            altered['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(altered)

    def test_omitted_obligation_or_failed_identity_rejected(self):
        altered = deepcopy(self.record)
        altered['remaining_obligations'].pop()
        with self.assertRaises(ValueError):
            validate(altered)
        altered = deepcopy(self.record)
        altered['reviews']['original_core_equations_and_forward_stress_jets']['rows'][0]['identities']['axial'] = False
        with self.assertRaises(ValueError):
            compare_records(self.record, altered)

    def test_cli_does_not_overwrite_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)/'record.json'
            target.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output', str(target)])
            self.assertEqual(target.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
