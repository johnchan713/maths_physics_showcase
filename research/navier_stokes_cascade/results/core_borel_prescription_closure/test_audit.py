#!/usr/bin/env python3
"""Check arbitrary-index algebra, rule identity, counterchecks and scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review
from proof import HERE,poly,add,certificates,verify_certificate,source_rule_checks,real_stage_residuals,floor_inputs_from_rule


class PrescriptionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        rows=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(rows);self.assertTrue(all(rows))

    def test_universal_normalized_margin_clears_positive_denominator(self):self.control('normalized_universal_polynomial')
    def test_universal_separation_and_boundary_margins(self):
        self.control('separation_universal_polynomial');self.control('boundary_universal_polynomial')

    def test_inherited_physical_log_power_row_is_preserved(self):self.control('physical_log_power_universal_polynomial')
    def test_wrong_symbolic_targets_are_rejected(self):self.control('wrong_target_is_rejected')
    def test_negative_remainder_cannot_certify_a_lower_margin(self):
        c=certificates()['normalized'];bad=dict(c,remainder=poly((-1,('g',))))
        bad['target']=add(bad['premise'],bad['remainder'])
        self.assertFalse(verify_certificate(bad))

    def test_polynomial_canonicalization_does_not_depend_on_variable_order(self):
        self.assertEqual(poly((1,('g','b')),(-1,('b','g'))),{})

    def test_historical_rule_contains_every_derivative_row_through_n(self):self.control('all_normalized_rows_through_n_required')
    def test_historical_rule_gain_is_n_h(self):self.control('fixed_gain_is_n_times_h')
    def test_historical_rule_keeps_separation_boundary_and_padding(self):self.control('original_max_includes_each_row')
    def test_dropped_derivative_row_and_wrong_gain_mutations_fail(self):
        self.control('missing_all_order_row_guard_changes');self.control('wrong_gain_changes')

    def test_missing_index_padding_boundary_and_separation_mutations_fail(self):
        for fragment in ('dropped_index_budget','dropped_unit_padding','dropped_domain_margin','dropped_separation'):self.control(fragment)

    def test_every_manufactured_stage_row_has_strict_gain_margin(self):self.control('every_finite_row_has_the_strict_gain_margin')
    def test_finite_prefix_is_not_an_infinite_certificate(self):
        self.control('no_actual_scale_selected_or_infinite_prefix_promotion');self.control('finite_prefix_does_not_certify_arbitrary_future')

    def test_later_conservative_bound_cannot_relabel_old_scale(self):self.control('larger_later_bound_is_not_automatically_satisfied')
    def test_invalid_or_float_log_standins_rejected(self):
        for args in ((0,1,1,0,0,0,1),(1,0,1,0,0,0,1),(1,1,0,0,0,0,1),(1,1,1,-1,0,0,1)):
            with self.assertRaises(ValueError):real_stage_residuals(*args)
        with self.assertRaises(TypeError):real_stage_residuals(1.0,1,1,0,0,0,1)

    def test_actual_normalized_constant_formula_is_unchanged_for_all_tested_norms(self):self.control('same_complete_tuple_constant')
    def test_Euler_conversion_and_extra_curl_derivative_remain(self):
        self.control('second_Euler_order_needs_Bell_conversion');self.control('extra_cutoff_derivative')

    def test_floor_receives_exact_same_Ctan_and_complete_transport(self):self.control('same_Ctan_and_complete_transport')
    def test_rule_defined_diagonal_and_diagnostic_identity_inputs_are_removed(self):
        self.control('independent_diagonal_is_discharged_only_for_rule_defined');self.control('diagnostic_cutoff_identity_is_not_needed')

    def test_summed_low_order_rate_still_qh_not_a_fabricated_q2h_bound(self):self.control('same_head_free_first_order_sum')
    def test_actual_independent_sequence_and_band_cannot_be_promoted(self):self.control('actual_independent_numeric_sequence_and_band_are_not_certified')
    def test_no_new_analytical_assumptions_and_four_records_bound(self):
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertEqual(len(self.record['dependency_evidence_sha256']),4)
        self.assertTrue(all(self.record['gates'].values()))

    def test_frozen_evidence_reproduces_exactly(self):self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)
    def test_historical_scope_and_remaining_obligations_cannot_be_changed(self):
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
