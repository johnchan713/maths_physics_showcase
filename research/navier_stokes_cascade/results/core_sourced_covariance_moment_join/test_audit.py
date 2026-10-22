#!/usr/bin/env python3
"""Generated covariance, full moments and weighted-domain regression checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,fixture
from budgets import P,through,moment_budget,generated_moments,next_pole_ledger,operator_coverage,orders


class GeneratedCovarianceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        rows=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(rows);self.assertTrue(all(rows))

    def test_all_tensor_entries_equal_difference_of_full_covariances(self):self.control('complete_tensor_matches_difference_of_two_covariances')
    def test_axial_and_all_six_entries_change(self):self.control('all_six_entries_change')
    def test_second_harmonic_square_survives_without_old_second_mode(self):self.control('absent_old_second_mode_still_has_nonzero_square')
    def test_l1_tensor_bound_keeps_cross_and_square(self):self.control('full_tensor_each_entry_has_l1_cross_plus_square_bound')
    def test_independent_complete_graph_curl_majorants(self):self.control('full_graph_cylindrical_curl_norm_is_bounded')
    def test_cutoff_curl_and_connections_cannot_be_erased(self):
        self.control('cutoff_derivatives_cannot_be_erased');self.control('graph_and_cylindrical_connection_are_nonzero')

    def test_full_radial_and_axial_covariance_moment_identities(self):
        self.control('full_axis_pressure_defect_identity');self.control('full_axial_covariance_and_pressure_moment')

    def test_compact_pressure_both_endpoints_and_full_first_moment(self):
        self.control('compact_pressure_both_endpoints');self.control('full_pressure_first_moment')

    def test_whole_moving_compensator_product_is_differentiated(self):
        self.control('whole_moving_pressure_derivative');self.control('moving_compensator_derivative_survives')

    def test_two_selected_fluxes_are_not_a_full_pressure_bound(self):self.control('axial_pressure_not_supplied_by_two_selected_fluxes')
    def test_independent_polynomial_moment_majorants(self):self.control('independent_polynomial_majorant_is_bounded')
    def test_generated_signed_modes_and_complete_cut_curl(self):self.control('all_signed_modes_and_complete_cut_curl_retained')
    def test_generated_covariance_is_not_an_independent_input(self):self.control('generated_covariance_cross_and_square_retained')
    def test_full_cut_ODE_pressure_is_retained(self):self.control('full_cut_pressure_keeps_unprojected_source')
    def test_all_generated_moments_are_supplied_by_one_way_chain(self):self.control('all_three_generated_moments_and_whole_pressure_response_are_positive')
    def test_finite_jet_dependency_reaches_before_initial_order_18(self):self.control('residual_order_zero_uses_finite_before_prefix_jets_through_18')
    def test_normal_bound_and_derived_tangent_are_used(self):self.control('normal_denominator_and_derived_tangent_still_used')
    def test_power_margins_keep_extra_axial_epsilon(self):self.control('all_power_margins_and_extra_axial_epsilon_retained')
    def test_new_pole_count_covers_self_square_and_fresh_postfast_rows(self):self.control('sourced_self_square_and_fresh_map_poles')
    def test_recomputed_operator_domains_pass(self):self.control('correctly_recomputed_domains_cover_next_maps')
    def test_old_operator_domains_rejected_for_positive_pole(self):self.control('old_domains_are_sufficient_only_at_zero_pole')
    def test_old_weighted_input_norm_does_not_bound_new_pole(self):self.control('old_4p_input_norm_cannot_bound_generic_6p_input')
    def test_same_width_primitive_recipe_exists_at_required_poles(self):self.control('same_fixed_width_primitive_recipe_exists_at_both_required_poles')

    def test_generated_moments_cannot_substitute_a_reserved_subhull(self):
        w=fixture()
        with self.assertRaises(TypeError):generated_moments(w,0,1,Q(1,2),lambda n:P(1),a=Q(3,2),d=2)

    def test_invalid_poles_hulls_dimensions_and_norms_rejected(self):
        for p in (-1,0.5):
            with self.assertRaises((ValueError,TypeError)):next_pole_ledger(p)
        for a,d in ((0,2),(2,1),(1,1)):
            with self.assertRaises(ValueError):moment_budget(0,1,1,a,d,1)
        with self.assertRaises(ValueError):moment_budget(0,-1,1,1,2,1)
        with self.assertRaises(ValueError):through(lambda n:P(1),0,0)
        with self.assertRaises(ValueError):operator_coverage(1,{})
        with self.assertRaises(ValueError):orders(0)

    def test_covariance_and_residual_orders_are_distinct(self):
        o=orders(Q(1,100000));self.assertEqual(o['covariance'],Q(3,2)-Q(1,100000))
        self.assertFalse(o['complete_residual_improvement_certified']);self.control('no_actual_enclosure_or_complete_residual_promotion')

    def test_four_frozen_records_and_no_added_assumptions(self):
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertEqual(self.record['protocol']['new_selector_or_smallness_conditions'],[])
        self.assertEqual(len(self.record['dependency_evidence_sha256']),4);self.assertTrue(all(self.record['gates'].values()))

    def test_new_frozen_evidence_reproduces_exactly(self):self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)
    def test_scope_or_obligations_cannot_be_promoted(self):
        for k,v in self.record['flags'].items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)

    def test_frozen_output_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
