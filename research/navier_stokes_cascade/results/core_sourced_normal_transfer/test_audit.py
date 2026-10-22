#!/usr/bin/env python3
"""Original floor transfer, exact denominator powers and pressure regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,fixture
from budgets import normal,phase,P,forced_from_original,degree_bound,transfer_factor


class SourcedNormalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        rows=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(rows);self.assertTrue(all(rows))

    def test_original_bound_enters_unchanged_inverse(self):self.control('same_original_normal_bound_enters_inverse')
    def test_stronger_phase_beta_is_not_inferred(self):self.control('stronger_optional_phase_beta_bound_is_not_inferred')
    def test_original_carriers_do_not_require_comparing_two_loose_bounds(self):self.control('original_carrier_bounds_without_artificial_cross_gate')
    def test_same_shell_tilt_growth_and_length_floor(self):self.control('original_shell_tilt_growth_and_L_ge_1_retained')
    def test_full_background_and_path_callbacks_are_retained(self):self.control('background_and_path_callbacks_are_same_objects')
    def test_source_tangent_and_length_inputs_are_retained(self):self.control('unprojected_source_and_tangent_length_are_same_objects')
    def test_all_full_phase_callbacks_reach_inverse(self):self.control('unchanged_full_phase_callbacks_reach_inverse')
    def test_no_actual_field_or_propagator_scope_promotion(self):self.control('no_actual_input_or_full_propagator_promotion')
    def test_formal_positive_q_degrees(self):self.control('positive_formal_q_degree')
    def test_nondegenerate_outputs_attain_degree_bound(self):self.control('top_degree_is_present')
    def test_formal_q_substitutions_match_original_recurrences(self):self.control('formal_matches_unchanged_inverse')
    def test_coefficientwise_scale_transfer_including_smaller_q(self):self.control('coefficientwise_scale_transfer')
    def test_scale_powers_are_unchanged(self):self.control('scale_powers_unchanged')
    def test_pure_normal_source_has_nonzero_pressure_when_velocity_is_zero(self):self.control('pure_normal_source_nonzero_pressure_with_zero_velocity')
    def test_moving_normal_survives(self):self.control('moving_normal_survives_pressure_and_ambient')
    def test_same_cutoff_pressure_and_potential_derivatives_survive(self):self.control('same_cutoff_derivative_products_survive')

    def test_degree_induction_and_repeated_ambient_cost_at_large_orders(self):
        for n in (1,2,10,1000):
            self.assertEqual(degree_bound('velocity',n),4*n+2)
            for i in (1,n):
                self.assertLessEqual(2+(i+1)+degree_bound('velocity',n-i),degree_bound('velocity',n))
            self.assertEqual(degree_bound('pressure',n),4*n+3)
        self.assertEqual(transfer_factor('velocity',2,2,8),4**10)

    def test_verifier_does_not_reselect_band_or_call_stronger_phase_gate(self):
        p,sp,ell,kw=fixture();old_select=normal.selector.select_tail;old_gate=phase.normal_gate
        def forbidden(*a,**k):raise AssertionError('Band reselection or stronger phase gate is forbidden')
        normal.selector.select_tail=forbidden;phase.normal_gate=forbidden
        try:self.assertTrue(forced_from_original(p,sp,ell,**kw)[2]['conditional_normal_floor_implication_passes'])
        finally:normal.selector.select_tail=old_select;phase.normal_gate=old_gate

    def test_early_band_shortened_length_and_later_enlarged_norm_rejected(self):
        p,sp,ell,kw=fixture()
        for pp,ss,ee in ((p,sp,16),(p,dict(sp,length=1),ell),(dict(p,Mz=10**30),sp,ell)):
            with self.assertRaises(ValueError):forced_from_original(pp,ss,ee,**kw)

    def test_missing_original_inputs_and_invalid_harmonic_rejected(self):
        p,sp,ell,kw=fixture();p.pop('Mz')
        with self.assertRaises(ValueError):forced_from_original(p,sp,ell,**kw)
        p,sp,ell,kw=fixture()
        with self.assertRaises(ValueError):forced_from_original(p,sp,ell,harmonic=0,**kw)

    def test_invalid_degree_and_scalar_transfer_inputs_rejected(self):
        for n in (-1,1.5,True):
            with self.assertRaises((ValueError,TypeError)):degree_bound('velocity',n)
        for q in (0,-1):
            with self.assertRaises(ValueError):transfer_factor('pressure',1,1,q)
        with self.assertRaises(ValueError):degree_bound('unknown',0)

    def test_three_historical_records_and_no_new_analytical_assumptions(self):
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertEqual(len(self.record['dependency_evidence_sha256']),3);self.assertTrue(all(self.record['gates'].values()))

    def test_new_frozen_record_reproduces_exactly(self):self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)
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
