#!/usr/bin/env python3
"""Axial determinant, full energy, scalar root and evidence integrity checks."""
from copy import deepcopy
from fractions import Fraction as F
import json
import tempfile
from pathlib import Path
import unittest
import mpmath as mp
from axial_bounds import exp_bounds,dyadic,matrix_review,shape_review,budget_review,prefix_tail_review
from axial_checks import exact_affine_review,numerical_matrix_review,scalar_review
from audit import HERE,FLAGS,generate,validate,compare_records,protocol,parent_review,main

class AxialRootTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_live_parent_and_protocol_match(self):self.assertEqual(protocol()['parent_commit'],'a41b941add8dc9c0368e3cea7c4f0aed979a1f78');self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
    def test_rational_scaled_inverse_rows_below_four(self):self.assertTrue(all(F(x)<4 for x in matrix_review()['scaled_inverse_row_upper']))
    def test_gap_lower_bound_keeps_positive_lambda(self):self.assertIn('2lambda',matrix_review()['exact_gap']);self.assertFalse(FLAGS['matrix_invertibility_at_lambda_zero_claimed'])
    def test_matrix_parameter_zero_is_singular_in_exact_fixture(self):
        from axial_checks import inverse_helpers
        with self.assertRaises(ValueError):inverse_helpers.inverse([[F(3,10),F(9,10)],[F(3,10),F(9,10)]])
    def test_trial_domain_includes_zero_for_constant_term_bound(self):self.assertIn('0<=A<=6/5',matrix_review()['target_C1'])
    def test_fixed_trial_and_composed_derivatives_differ(self):self.assertTrue(exact_affine_review()['checks']['composed_coefficient_derivative_requires_strength_chain_term'])
    def test_positive_Taylor_proves_binary_majorant(self):self.assertGreater(exp_bounds(8)[0],2**11);self.assertLess(F(64000)*F(2)**-352,F(1,10**100))
    def test_dyadic_conversion_is_exact_for_both_signs(self):
        self.assertEqual(dyadic(mp.mpf('.125')),F(1,8));self.assertEqual(dyadic(mp.mpf('-1.25')),F(-5,4));self.assertEqual(dyadic(mp.mpf(0)),0)
    def test_monotone_shape_enclosure_has_1024_boxes(self):self.assertEqual(shape_review()['boxes'],1024);self.assertEqual(shape_review()['interval_digits'],80)
    def test_shape_rationals_are_inside_required_positive_K_range(self):
        a,b=map(F,shape_review()['rational_K_enclosure']);self.assertTrue(F(49,200)<a<b<F(1,4))
    def test_independent_first_ramp_bound_is_sufficient(self):self.assertLess(241+F(23,6*64),500);self.assertIn('11/6',prefix_tail_review()['first_energy_mass_bound'])
    def test_prefix_bounds_keep_full_eta_norm(self):self.assertIn('||m0||_C1',prefix_tail_review()['prefix']);self.assertLess(24*3*F(2)**-7712,1)
    def test_complete_post_mass_fits_256T_at_smallest_T(self):self.assertLess(F(1002,64)+240,256);self.assertIn('1024T',prefix_tail_review()['full_post_mass'])
    def test_release_bound_retains_infinite_tail(self):self.assertLess(1+F(1,2)+F(2,3)*F(1,100)**7,2);self.assertIn('nonzero end bump',prefix_tail_review()['release_mass'])
    def test_eta_density_includes_angular_edit_derivative(self):self.assertLess(2+2*F(1,100)/F(99,100),4)
    def test_superexponential_absorption_is_logarithmic(self):self.assertGreater(exp_bounds(7)[0],800);self.assertLess((7+21*256)*F(2)**-352,F(4,5))
    def test_left_and_right_root_signs_recomputed_with_rationals(self):
        a,b=map(F,shape_review()['rational_K_enclosure']);el,eh=exp_bounds(-26,128);err=64000*F(2)**-352
        self.assertLess(F('1.0100502')**2*b-(1-eh)/4+err,F(-3,10**8))
        self.assertGreater(F('1.0100504')**2*a-(1-el)/4-err,F(6,10**8))
    def test_same_strength_bracket_is_preserved(self):self.assertEqual(budget_review()['amplitude_bracket'],['1.0100502','1.0100504'])
    def test_positive_discriminant_proves_root_jacobian(self):self.assertGreater(4*F(49,200)*F(249,1000),F(49,100)**2);self.assertEqual(budget_review()['root_inverse_upper'],'100/49')
    def test_sharper_eta_bound_follows_from_complete_error(self):
        ratio=F(1,64)+512*exp_bounds(-26,128)[1]+80*F(2)**(-352*40)/64
        self.assertLess(ratio,F(49,3000));self.assertEqual(budget_review()['root_eta_bound'],'T*lambda/30')
    def test_end_energy_has_nonzero_positive_increments(self):self.assertTrue(exact_affine_review()['checks']['positive_quadratic_and_constant_increments_retained'])
    def test_full_affine_energy_cross_term_is_not_erased(self):self.assertTrue(exact_affine_review()['checks']['full_affine_end_energy_identity_exact']);self.assertTrue(exact_affine_review()['checks']['cross_increment_nonzero'])
    def test_odd_cross_term_prevents_assumed_even_strength(self):self.assertTrue(exact_affine_review()['checks']['strength_not_assumed_even_in_eta']);self.assertFalse(FLAGS['strength_evenness_assumed'])
    def test_both_normalized_moments_close_exactly(self):self.assertEqual(exact_affine_review()['exact_proxy_moment_components'],4)
    def test_correcting_only_M_leaves_J_nonzero(self):self.assertTrue(exact_affine_review()['checks']['only_M_correction_leaves_nonzero_J'])
    def test_three_actual_bump_matrix_diagnostics(self):
        rows=numerical_matrix_review()['rows'];self.assertEqual(len(rows),3);self.assertTrue(all(r['scaled_inverse']<4 and r['absolute_row_error']<1e-60 for r in rows))
    def test_diagnostic_parameters_do_not_evaluate_selected_lambda(self):self.assertTrue(numerical_matrix_review()['finite_quadratures_are_diagnostics']);self.assertTrue(numerical_matrix_review()['selected_lambda_not_evaluated'])
    def test_actual_factored_solver_beats_collapsed_binary_gap(self):self.assertTrue(numerical_matrix_review()['checks']['binary64_subtraction_can_erase_nonzero_gap']);self.assertTrue(numerical_matrix_review()['checks']['actual_factored_solver_matches_independent_LU'])
    def test_scalar_root_and_implicit_derivative_controls(self):
        r=scalar_review();self.assertTrue(r['checks']['exact_rational_discriminant_jacobian_identity']);self.assertTrue(all(v['derivative_error']<1e-70 for v in r['rows']))
    def test_complete_heat_net_energy_identity_is_retained(self):self.assertTrue(FLAGS['heat_and_compensation_exact_net_S_identity_retained']);self.assertIn('renormalized-I',protocol()['heat'])
    def test_all_fixed_jets_do_not_remove_other_construction_dependencies(self):self.assertTrue(FLAGS['all_fixed_one_sided_root_and_coefficient_jets_supplied']);self.assertFalse(FLAGS['all_retained_analytic_dependencies_discharged'])
    def test_all_1051_historical_and_18_direct_bindings(self):self.assertEqual(len(self.record['inherited_provenance']),1051);self.assertEqual(len(self.record['source_hashes']),18)
    def test_frozen_record_and_json_roundtrip_match(self):
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.record,f);validate(json.loads(json.dumps(self.record)))
    def test_unsupported_scope_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_history_key_substitution_and_cached_mutation_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        r=generate();r['reviews']['actual_scaled_axial_inverse']['checks']['zero_lambda_not_invertible']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1]);self.assertTrue(matrix_review()['checks']['zero_lambda_not_invertible'])
    def test_protocol_sources_gates_scope_and_overwrite_guards(self):
        for k in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[k],dict):r[k].pop(next(iter(r[k])))
            else:r[k].pop()
            with self.assertRaises(ValueError):validate(r)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')

if __name__=='__main__':unittest.main()
