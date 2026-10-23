#!/usr/bin/env python3
"""Current scalar-gain checks; no finite test is an actual-scale enclosure."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review
from budgets import P,rates,slow_reset,fast_change


class ActualScalarGainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        a=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(a);self.assertTrue(all(a))

    def test_independent_unreduced_plain_pressure_budget(self):self.control('independent_general_dynamic_P')
    def test_independent_weighted_radial_source_budget(self):self.control('independent_general_dynamic_Qrem')
    def test_independent_angular_scalar_budget(self):self.control('independent_general_dynamic_Jtheta')
    def test_independent_axial_scalar_budget(self):self.control('independent_general_dynamic_Jz')
    def test_subunit_B_plus_one_error_detected(self):self.control('subunit_square_loses_B_plus_one')
    def test_affine_rate_bound_covers_whole_interval(self):self.control('whole_interval_')
    def test_strict_C_margin_on_original_kappa_range(self):self.control('coarse_scalar_margin_strict')
    def test_scalar_gain_does_not_require_fast_D(self):self.control('current_scalars_reach_C_without_fast_D')
    def test_fast_feedback_still_blocks_full_gain(self):self.control('fast_feedback_still_blocks_D')
    def test_nonzero_G_hull_leak_in_angular_scalar(self):self.control('support_contract_violation_leaks_linear_angular_scalar')
    def test_nonzero_G_hull_leak_in_axial_scalar(self):self.control('support_contract_violation_leaks_linear_axial_scalar')
    def test_leak_witness_keeps_both_velocity_moments(self):self.control('witness_both_reserved_velocity_moments_zero')
    def test_whole_moving_pressure_product_required(self):
        self.control('whole_moving_pressure_product_derivative');self.control('frozen_compensator_derivative_is_detectably_wrong')

    def test_positive_swirl_square_is_not_higher_order(self):self.control('compact_swirl_square_pressure_strictly_positive')
    def test_asymmetric_complete_fast_components(self):self.control('independent_asymmetric_fast_')
    def test_every_old_wave_and_fast_scalar_contribution_survives(self):
        for fragment in ('actual_old_reset_wave_fast_all_retained_','omitting_old_reset_changes_',
                         'omitting_complete_covariance_changes_','omitting_actual_fast_correlations_changes_'):self.control(fragment)

    def test_fresh_slow_uses_extra_current_angular_jet(self):self.control('fresh_slow_uses_derived_current_scalars_and_extra_angular_jet')
    def test_whole_axial_pressure_moment_composes(self):self.control('whole_axial_moment_keeps_full_current_pressure')
    def test_same_pole_recipes_and_long_shell_embeddings(self):
        self.control('same_reweighted_original_operator_recipes');self.control('long_shell_lower_pole_embedding_is_retained')

    def test_A1_and_original_enclosure_limits_retained(self):
        self.control('same_complete_fast_A1_field_budget');self.control('original_actual_input_enclosures_not_promoted')

    def test_full_source_gain_interface_is_not_opened(self):self.control('improved_full_source_interface_still_blocked')
    def test_finite_dependency_ledger(self):self.assertEqual(self.r['parts']['composition']['highest_original_initial_jet'],26)
    def test_invalid_ranges_and_missing_components(self):
        for k in (0,Q(1,4),-1,0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
        for fn,args in ((slow_reset,({}, {},P(1),P(1),P(1))),(fast_change,({}, {},P(1)))):
            with self.assertRaises(ValueError):fn(0,*args)

    def test_frozen_record_and_seven_dependencies(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),7)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))

    def test_scope_and_obligation_mutation_rejected(self):
        for k,v in self.record['flags'].items():
            r=deepcopy(self.record);r['flags'][k]=not v
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)

    def test_historical_hash_mutation_rejected(self):
        r=deepcopy(self.record);r['dependency_evidence_sha256'][next(iter(r['dependency_evidence_sha256']))]='bad'
        with self.assertRaises(ValueError):audit.validate(r)

    def test_existing_frozen_file_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
