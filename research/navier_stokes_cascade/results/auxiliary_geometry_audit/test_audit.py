#!/usr/bin/env python3
"""All-band geometry, exact covering certificates, source checks and scope."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from geometry import (COLORS,PRIME,J,DELTA_MAX,adjugate,center,color,coloring_review,counting_review,
                       covering_gap_review,covering_review,determinant,kernel_points,matrix_power,
                       multiply,norm_infinity,prime_certificate,radius_bound,separation_review)
from checks import (KAPPA,band_logs,chain_rule_review,compatibility_review,failed_shortcuts_review,
                     partition_review,path_review,scale_and_source_review,scale_exponent_certificate)
from audit import FLAGS,HERE,compare_records,expected_protocol,inherited_provenance,main,validate


class ExactGeometryTests(unittest.TestCase):
    def test_correct_integer_matrix_and_covering_powers(self):
        self.assertEqual(J,((3,1),(1,5)))
        self.assertEqual(matrix_power(2),((10,8),(8,26)))
        self.assertEqual(matrix_power(3),((38,50),(50,138)))
        self.assertEqual(norm_infinity(matrix_power(3)),188)

    def test_adjugates_are_exact_for_every_common_covering(self):
        for d in range(DELTA_MAX+1):
            a=matrix_power(d)
            det=determinant(a)
            self.assertEqual(multiply(a,adjugate(a)),((det,0),(0,det)))

    def test_all_band_covering_difference_bound_uses_rational_majorants(self):
        row=covering_gap_review()
        self.assertEqual(row['rational_covering_difference_upper_bound'],F(1873,625))
        self.assertTrue(all(row['checks'].values()))

    def test_complete_residue_coloring_has_1250_values(self):
        row=coloring_review()
        self.assertEqual(row['count'],1250)
        self.assertTrue(all(row['checks'].values()))

    def test_negative_grid_indices_and_translates_obey_the_same_color_rule(self):
        self.assertEqual(color(16,(-2,7,-11),1),color(21,(3,-3,4),1))
        self.assertNotEqual(color(16,(-2,7,-11),0),color(16,(-2,7,-11),1))

    def test_prime_and_all_relevant_center_covering_constraints(self):
        self.assertTrue(prime_certificate(PRIME))
        self.assertGreater(PRIME,COLORS)
        self.assertEqual(center(1),(0,F(1,2003)))
        row=separation_review()
        self.assertEqual(row['off_diagonal_entries'],[1,8,50])
        self.assertTrue(all(row['checks'].values()))

    def test_exact_closed_rectangle_radius_and_congruence_margin(self):
        self.assertEqual(radius_bound(),F(1,6057072))
        row=separation_review()
        self.assertEqual(row['maximum_enlarged_difference_displacement'],F(1,8012))
        self.assertEqual(row['center_congruence_separation_margin'],F(3,8012))

    def test_common_torus_kernel_points_are_all_present(self):
        for delta,expected in enumerate((1,14,196,2744)):
            denominator,points=kernel_points(delta)
            self.assertEqual(denominator,expected)
            self.assertEqual(len(points),expected)
        self.assertTrue(all(covering_review()['checks'].values()))

    def test_pointwise_overlap_and_radial_counts_are_polynomial(self):
        row=counting_review()
        self.assertEqual(row['original_pointwise_slow_label_bound'],162)
        self.assertEqual(row['common_lift_factor'],2744)
        self.assertTrue(all(row['checks'].values()))

    def test_invalid_band_color_and_covering_parameters_are_rejected(self):
        for args in ((15,(0,0,0),0),(16,(0.,0,0),0),(16,(0,0,0),2)):
            with self.assertRaises(ValueError):
                color(*args)
        for n in (0,1251):
            with self.assertRaises(ValueError):
                center(n)
        with self.assertRaises(ValueError):
            kernel_points(4)


class SourceAndChartTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scales=scale_and_source_review()
        cls.partition=partition_review()
        cls.paths=path_review()
        cls.chain=chain_rule_review()

    def test_source_radial_exponent_has_outer_parentheses(self):
        self.assertEqual(self.scales['canonical_rendered_formula'],'d_r=2*((1+h)*rho_g-h*kappa_s)')
        self.assertTrue(self.scales['checks']['frozen_fixture_agrees_with_visually_verified_formula'])

    def test_suspected_source_error_is_explicitly_rejected_not_promoted(self):
        self.assertTrue(FLAGS['suspected_factor_two_normalization_error_rejected'])
        self.assertTrue(self.scales['checks']['suspected_factor_two_error_is_explicitly_rejected'])
        self.assertLess(self.scales['rejected_variant_normalized_log_ratio'],-1000)

    def test_exact_epsilon_power_cancellation_is_not_a_band_scan(self):
        for h in (F(1,200),F(9,1000),F(1,1000000)):
            c=scale_exponent_certificate(h)
            self.assertEqual(c['cancelled_rho_coefficient'],0)
            self.assertEqual(c['remaining_epsilon_exponent'],-KAPPA)
            self.assertEqual(c['misparsed_epsilon_exponent'],-KAPPA/2)

    def test_huge_bands_match_floor_intervals_without_materializing_Q(self):
        self.assertEqual(max(r['ell'] for r in self.scales['rows']),10**12)
        self.assertTrue(self.scales['checks']['all_sampled_time_scales_match_exact_floor_interval'])
        self.assertTrue(self.scales['checks']['all_sampled_radial_scales_match_exact_floor_interval'])

    def test_partition_denominator_and_squared_coverage(self):
        self.assertTrue(all(self.partition['checks'].values()))

    def test_partition_is_not_double_counted_for_two_pulse_signs(self):
        self.assertTrue(self.partition['checks']['squared_partition_is_counted_once_per_box_not_twice_per_sign'])

    def test_pulse_paths_keep_transverse_coordinates_and_each_copy(self):
        for key in ('path_keeps_same_transverse_coordinate','every_lifted_path_has_the_correct_band_point',
                    'pulse_coordinate_advances_exactly'):
            self.assertTrue(self.paths['checks'][key])

    def test_path_derivative_and_time_speed_are_uniform_in_covering_difference(self):
        self.assertTrue(self.paths['checks']['path_derivative_is_the_fixed_projector'])
        self.assertTrue(self.paths['checks']['common_torus_time_speed_is_one'])

    def test_actual_evaluated_map_has_all_normalized_physical_derivatives(self):
        self.assertTrue(self.chain['checks']['actual_physical_jets_match_all_three_normalized_operators'])

    def test_missing_fast_graph_terms_are_detected(self):
        self.assertTrue(self.chain['checks']['dropping_fast_radial_or_temporal_terms_fails'])

    def test_band_and_h_validation(self):
        for args in ((15,F(1,200)),(16,0),(16,F(1,100))):
            with self.assertRaises(ValueError):
                band_logs(*args)


class CompatibilityAndFailureTests(unittest.TestCase):
    def test_specified_absolute_operators_preserve_pullbacks(self):
        self.assertTrue(all(compatibility_review()['checks'].values()))

    def test_distinct_centers_do_not_suffice_under_different_coverings(self):
        row=failed_shortcuts_review()
        self.assertTrue(row['checks']['distinct_band_centers_alone_do_not_imply_absolute_preimage_separation'])

    def test_period_four_fails_for_closed_doubled_boxes(self):
        self.assertTrue(failed_shortcuts_review()['checks']['period_four_fails_the_chosen_closed_doubled_box_graph'])

    def test_excessive_radius_can_create_overlap(self):
        self.assertTrue(failed_shortcuts_review()['checks']['radius_must_be_tied_to_covering_distances'])

    def test_common_sources_are_not_averaged_across_deck_copies(self):
        self.assertTrue(failed_shortcuts_review()['checks']['common_sources_must_not_be_averaged_over_other_preimage_copies'])

    def test_congruence_gap_is_not_promoted_to_an_absolute_torus_distance(self):
        self.assertTrue(failed_shortcuts_review()['checks']['congruence_margin_is_not_claimed_as_uniform_absolute_torus_Euclidean_distance'])


class EvidenceTests(unittest.TestCase):
    def test_protocol_and_all_inherited_source_bytes(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),expected_protocol())
        checks=inherited_provenance()
        self.assertEqual(len(checks),404)
        self.assertTrue(all(checks.values()))

    def test_geometric_proof_does_not_certify_analytic_foundation_or_blowup(self):
        row=json.loads((HERE/'evidence.json').read_text())
        validate(row)
        self.assertEqual(len(row['gates']),13)
        self.assertFalse(row['flags']['full_global_pulse_and_nonlinear_correction_hierarchy_revalidated'])
        self.assertFalse(row['flags']['blowup_verified'])

    def test_rejected_false_alarm_cannot_be_silently_reintroduced(self):
        row=json.loads((HERE/'evidence.json').read_text())
        row['flags']['suspected_factor_two_normalization_error_rejected']=False
        with self.assertRaises(ValueError):
            validate(row)

    def test_unconditional_promotion_is_rejected(self):
        for flag in ('blowup_verified','unconditional_blowup_construction_verified',
                     'leading_profile_background_and_all_order_inverse_foundation_independently_verified'):
            row=json.loads((HERE/'evidence.json').read_text())
            row['flags'][flag]=True
            with self.assertRaises(ValueError):
                validate(row)

    def test_exact_radius_or_compatibility_scope_tampering_is_rejected(self):
        row=json.loads((HERE/'evidence.json').read_text())
        bad=copy.deepcopy(row)
        bad['reviews']['rational_centers_and_closed_separation']['radius_upper_bound']='1/2003'
        with self.assertRaises(ValueError):
            compare_records(row,bad)
        bad=copy.deepcopy(row)
        bad['geometric_proof']['analytic_leading_profile_and_global_PDE_hierarchy_not_reproved']=False
        with self.assertRaises(ValueError):
            validate(bad)

    def test_frozen_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'replay.json'
            main(['--verify-record',str(HERE/'evidence.json'),'--output',str(output)])
            self.assertTrue(output.exists())
            with self.assertRaises(SystemExit):
                main(['--output',str(output)])


if __name__=='__main__':
    unittest.main()
