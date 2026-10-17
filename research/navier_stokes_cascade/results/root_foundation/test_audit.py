#!/usr/bin/env python3
"""Actual angular/heat maps, full C1 contraction, terminal wait and scope."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from root_bounds import exp_bounds,exp_interval,outward,mass_interval,matrix_review,binary_upper,root_budget,derivative_bounds,target_review,smallness_review,waiting_review
from root_checks import inverse,lagrange_inverse,product,identity,exact_map_review,numerical_bump_review,jets_and_wait_review
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main

class RootFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_parent_is_exact(self):self.assertEqual(protocol()['parent_commit'],'bf4794ad6e0dc309acab0cc6e84ad74a8ea558a1')
    def test_positive_Taylor_exponential_enclosure(self):
        a,b=exp_bounds(1);self.assertLess(a,b);self.assertGreater(a,F(27,10));self.assertLess(b,F(11,4))
    def test_negative_exponential_is_reciprocal_enclosure(self):
        a,b=exp_bounds(2);self.assertEqual(exp_bounds(-2),(1/b,1/a))
    def test_zero_exponential_is_exact(self):self.assertEqual(exp_bounds(0),(1,1))
    def test_exponential_invalid_order_or_range_rejected(self):
        for x,n in ((1,0),(1,True),(66,64)):
            with self.assertRaises(ValueError):exp_bounds(x,n)
        with self.assertRaises(ValueError):exp_interval(1,0)
    def test_outward_rounding_includes_negative_and_positive_rationals(self):
        for q in (F(1,3),F(-2,7),F(3)):
            a,b=outward(q,100);self.assertLessEqual(a,q);self.assertGreaterEqual(b,q)
    def test_bump_mass_uses_full_width_and_pressure_factor_two(self):
        self.assertEqual(mass_interval((0,0)),(F(3,10),F(3,10)));self.assertEqual(mass_interval((0,0),factor=2),(F(3,5),F(3,5)))
    def test_bad_bump_data_rejected(self):
        for args in (((1,0),0,1,1),((0,1),0,0,1),((0,1),0,1,0)):
            with self.assertRaises(ValueError):mass_interval(*args)
    def test_uniform_angular_inverse_below_three(self):self.assertTrue(all(F(x)<3 for x in matrix_review()['angular_inverse_row_upper']))
    def test_uniform_heat_inverse_below_fifteen(self):self.assertTrue(all(F(x)<15 for x in matrix_review()['heat_inverse_row_upper']))
    def test_uniform_quadratic_norms_include_square_terms(self):self.assertLess(F(matrix_review()['angular_quadratic_upper']),32);self.assertEqual(F(matrix_review()['heat_quadratic_upper']),F(144,5))
    def test_Lagrange_coefficients_match_independent_inverse(self):
        nodes=[F(1,8),F(1),F(7)];masses=[F(1,6),F(1,3),F(1,2)];B=[[m*n**j for j in range(3)] for n,m in zip(nodes,masses)]
        J=lagrange_inverse(nodes,masses);self.assertEqual(J,inverse(B));self.assertEqual(product(B,J),identity(3))
    def test_coincident_nodes_or_singular_matrices_rejected(self):
        with self.assertRaises(ValueError):lagrange_inverse([F(1),F(1)],[F(1),F(1)])
        with self.assertRaises(ValueError):inverse([[1,1],[1,1]])
    def test_44_exact_inverse_product_components_and_transpose_control(self):self.assertEqual(exact_map_review()['exact_two_sided_matrix_components'],44);self.assertTrue(exact_map_review()['checks']['transposed_inverse_detected'])
    def test_relative_square_moment_row_normalizations_exact(self):
        r=exact_map_review()['checks'];self.assertTrue(r['angular_relative_square_keeps_two_linear_plus_quadratic']);self.assertTrue(r['heat_pressure_and_minus_S_keep_linear_plus_half_quadratic']);self.assertTrue(r['I_row_is_exactly_linear'])
    def test_omitted_quadratic_energy_and_pressure_terms_nonzero(self):self.assertTrue(exact_map_review()['checks']['omitted_quadratic_terms_are_nonzero'])
    def test_six_actual_bump_integral_diagnostics(self):self.assertEqual(len(numerical_bump_review()['rows']),6);self.assertTrue(numerical_bump_review()['checks']['independent_integrals_match_exact_factoring'])
    def test_finite_quadrature_samples_not_continuum_bounds(self):self.assertTrue(numerical_bump_review()['finite_sample_integrals_not_continuum_enclosures']);self.assertTrue(numerical_bump_review()['actual_selected_lambda_and_roots_not_evaluated'])
    def test_matrix_range_not_promoted_to_target_range(self):self.assertIn('exp(-256)',target_review()['angular_target']);self.assertFalse(FLAGS['matrix_parameter_range_promoted_to_target_range'])
    def test_heat_normalization_keeps_both_inverse_f_factors(self):self.assertEqual(target_review()['heat_normalized_row_ratios_to_common_bound'],['27/80000000000','27/40000000000','1'])
    def test_finite_h_infinite_angular_integral_kept(self):self.assertTrue(target_review()['checks']['infinite_heat_angular_integral_h_times_one_over_h_retained'])
    def test_strict_full_C1_root_ball_and_half_contraction(self):
        r=root_budget(3,32,F(1,10**8));self.assertLess(r['lipschitz'],F(1,2));self.assertLessEqual(r['image_upper'],r['radius']);self.assertEqual(r['jacobian_inverse_upper'],6)
    def test_zero_target_and_zero_quadratic_cases(self):self.assertEqual(root_budget(3,32,0)['radius'],0);self.assertEqual(root_budget(3,0,1)['lipschitz'],0)
    def test_failed_contraction_is_rejected(self):
        with self.assertRaises(ValueError):root_budget(3,32,1)
        with self.assertRaises(ValueError):root_budget(0,1,1)
    def test_binary_majorants_do_not_materialize_selected_lambda(self):self.assertEqual(binary_upper(1,30),-7680);self.assertEqual(binary_upper(F(1,4),0),-2)
    def test_sharper_roots_fit_original_uniqueness_balls(self):self.assertEqual(smallness_review()['angular']['root_coefficient'],4200);self.assertEqual(smallness_review()['heat']['root_coefficient'],1920);self.assertTrue(smallness_review()['checks']['sharper_roots_lie_inside_original_uniqueness_balls'])
    def test_both_relative_edit_and_radial_slope_gaps(self):
        for k in ('angular','heat'):
            self.assertLess(smallness_review()[k]['relative_edit_log2_upper'],-7);self.assertLess(smallness_review()[k]['slope_over_lambda_log2_upper'],-1)
    def test_constant_Q_normalized_derivative_recursion(self):self.assertEqual(derivative_bounds(1,1,[0,1,0,0],0),[0,2,8,64])
    def test_twelve_Catalan_jet_controls_and_factorial_failure(self):self.assertEqual(jets_and_wait_review()['exact_scalar_series_orders'],12);self.assertTrue(jets_and_wait_review()['checks']['second_derivative_factorial_is_not_normalized_coefficient'])
    def test_heat_one_sided_jets_do_not_claim_complex_strip(self):self.assertFalse(FLAGS['common_complex_strip_for_heat_targets_claimed'])
    def test_terminal_wait_has_strict_input_and_target_bounds(self):self.assertGreater(F(waiting_review()['q_before_lower']),1);self.assertTrue(0<F(waiting_review()['qp_upper'])<1)
    def test_terminal_integrating_factor_and_nonzero_prescribed_constant(self):self.assertTrue(jets_and_wait_review()['checks']['cleared_terminal_integrating_factor_identity_exact']);self.assertTrue(jets_and_wait_review()['checks']['waiting_target_cannot_be_replaced_by_zero'])
    def test_waiting_is_scalar_eta_independent_and_strength_independent(self):self.assertTrue(waiting_review()['checks']['unique_wait_is_positive_finite_and_eta_independent']);self.assertTrue(waiting_review()['checks']['waiting_does_not_depend_on_the_axial_strength_root'])
    def test_all_1033_inherited_and_18_direct_bindings(self):self.assertEqual(len(self.record['inherited_provenance']),1033);self.assertEqual(len(self.record['source_hashes']),18)
    def test_frozen_record_and_roundtrip_match(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r);validate(json.loads(json.dumps(self.record)))
    def test_other_root_and_unconditional_scope_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_cached_reviews_not_mutable_aliases(self):
        r=generate();r['reviews']['exact_matrix_and_moment_map_controls']['checks']['I_row_is_exactly_linear']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(exact_map_review()['checks']['I_row_is_exactly_linear'])
    def test_same_count_historical_key_substitution_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])
    def test_protocol_gates_sources_and_other_obligations_cannot_be_erased(self):
        for k in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[k],dict):r[k].pop(next(iter(r[k])))
            else:r[k].pop()
            with self.assertRaises(ValueError):validate(r)
    def test_cli_refuses_evidence_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')

if __name__=='__main__':unittest.main()
