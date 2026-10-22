#!/usr/bin/env python3
"""Meaningful exact checks for rate/order limits and the full axial floor."""
from copy import deepcopy
from fractions import Fraction as Q
import inspect
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import diagonal_tail,comparison,tangential_error,normal_gate,length_compatibility
from controls import fixture,review


class TangentialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        rows=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(rows);self.assertTrue(all(rows))

    def test_partial_sum_and_infinite_remainder_identity(self):self.control('partial_plus_exact_remainder')
    def test_both_low_orders_have_universal_qh_majorant(self):self.control('both_low_orders_bounded_by_t')
    def test_low_order_unit_constant_is_sharp_for_majorant(self):self.control('unit_constant_is_sharp')
    def test_new_majorant_does_not_claim_q2h_without_head(self):self.control('cannot_give_uniform_t_squared')
    def test_higher_order_diagonal_gap_is_rejected(self):
        for m in (-1,2,3,True):
            with self.assertRaises(ValueError):diagonal_tail(Q(1,2),m)

    def test_invalid_geometric_series_inputs_rejected(self):
        for t in (0,-1,2):
            with self.assertRaises(ValueError):diagonal_tail(t)
        with self.assertRaises(ValueError):diagonal_tail(1,start=0)
        with self.assertRaises(TypeError):diagonal_tail(0.5)

    def test_all_saturated_slow_allocations_have_Ctan_over_S_bound(self):self.control('original_slow_and_component_sum')
    def test_complete_axial_transport_survives(self):self.control('full_axial_transport_is_not_dropped')
    def test_nearest_nonzero_rounding_uses_unit_cost(self):self.control('nonzero_rounding_cost_is_one_not_half')
    def test_tilt_and_box_geometry_survive(self):self.control('both_tilt_and_box_costs_survive')
    def test_original_epsilon_squared_gate_retained(self):self.control('epsilon_squared_remains')
    def test_failed_slow_gate_cannot_be_promoted(self):self.control('slow_failure_is_reported')
    def test_consistent_squared_epsilon_and_actual_ceiling_carrier_relation(self):self.control('consistent_epsilon_squared_and_ceiling_carrier')
    def test_old_full_recipe_is_strictly_larger_across_first_cutoff_budgets(self):self.control('old_full_Cphi_strictly_larger')
    def test_axial_transport_improvement_is_exact(self):self.control('old_axial_excess_exact')
    def test_new_margin_can_pass_when_old_sufficient_test_fails(self):self.control('same_supplied_scale_new_margin_passes_old_margin_fails')
    def test_new_margin_boundary_and_full_inverse_norm(self):self.control('new_exact_margin_retains_inverse_normal_square')
    def test_actual_parent_leading_recipe_requires_no_first_repair_for_this_join(self):self.control('existing_support_aware_leading_first_rows_feed')
    def test_original_length_factor_matches_direct_normal_gate(self):self.control('original_radius_and_length_constants_match')
    def test_insufficient_original_length_reported_without_reselection(self):
        self.control('insufficient_existing_length_is_reported');self.control('does_not_select_scale')

    def test_no_new_independent_low_order_inputs(self):
        self.control('no_independent_first_K2_or_leading_second_input')
        self.assertEqual(set(inspect.signature(comparison).parameters),{'Lead1','Rmin','Rmax','U','lambda_max','gmin','Llow','Lhigh','Dbox'})
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])

    def test_infinite_diagonal_premise_is_not_inferred_from_controls(self):self.control('original_infinite_diagonal_remains')
    def test_consistent_growth_and_nonnegative_transport_required(self):
        r=fixture()
        for lam in (0,4):
            with self.assertRaises(ValueError):normal_gate(r,lam,10,0,0,0)
        with self.assertRaises(ValueError):tangential_error(r,0,0,0,0)
        with self.assertRaises(ValueError):tangential_error(r,1,0,-1,0)
        with self.assertRaises(ValueError):length_compatibility(r,1,0,1)
        with self.assertRaises(ValueError):comparison(-1,1,2,2,3,1,1,2,6)
        with self.assertRaises(ValueError):comparison(1,2,1,2,3,1,1,2,6)

    def test_four_historical_records_and_declared_sources_match(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),4)
        self.assertTrue(all(self.record['gates'].values()))

    def test_frozen_record_reproduces_exactly(self):self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)
    def test_scope_and_remaining_obligations_cannot_be_promoted(self):
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
