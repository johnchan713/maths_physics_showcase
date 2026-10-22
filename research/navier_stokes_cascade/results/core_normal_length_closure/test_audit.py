#!/usr/bin/env python3
"""Original length closure, complete transport and stronger-bound distinction."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,fixture
from budgets import frame,selector,normal_length_certificate,tangential_error,verify_existing


class NormalLengthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        rows=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(rows);self.assertTrue(all(rows))

    def test_original_Riccati_budget_uniformly_dominates_normal_threshold(self):self.control('physical_consistency_and_uniform_length_implication')
    def test_CE_is_an_independent_primitive_formula_before_ODE_estimates(self):self.control('CE_is_primitive_positive_budget_before_any_ODE')
    def test_both_algebraic_factors_are_retained(self):self.control('both_strict_domination_factors_retained')
    def test_25_and_Riccati_parts_cover_all_normal_threshold_terms(self):self.control('original_normal_threshold_dominated')
    def test_reference_growth_identity_implies_consistency(self):self.control('actual_reference_identity_implies_consistency')
    def test_positive_but_unphysical_reference_data_are_not_certified(self):self.control('positive_unphysical_reference_inputs_are_not_promoted')
    def test_each_slow_budget_allocation_keeps_full_error_bound(self):self.control('full_error_bounded_by_original_Dt_over_L')
    def test_consistent_physical_epsilon_square_and_carrier_relation(self):self.control('consistent_epsilon_square_and_ceiling_carrier_relation')
    def test_full_axial_transport_uses_actual_L(self):self.control('complete_axial_transport_keeps_full_Hz_and_actual_L')
    def test_both_tilt_terms_and_radial_displacement_remain(self):self.control('both_tangential_tilt_terms_and_box_displacement_survive')
    def test_failed_slow_gate_is_reported(self):self.control('slow_failure_is_reported')
    def test_same_original_selector_closes_floor_for_all_supplied_later_bands(self):self.control('same_original_selector_and_length_close_normal_nonvanishing')
    def test_original_radius_and_existing_band_are_retained(self):self.control('actual_radius_factor_and_supplied_band_retained')
    def test_old_length_semantics_and_unknown_early_band_cannot_be_promoted(self):
        self.control('old_length_without_its_declared_semantics_is_rejected_not_reselected');self.control('unknown_early_band_is_not_identified')

    def test_new_upper_norm_cannot_be_silently_inserted_into_old_length(self):self.control('new_conservative_Mz_cannot_be_silently_substituted')
    def test_separately_valid_box_enclosures_need_no_extra_comparison_gate(self):self.control('independently_valid_box_enclosures_need_no_artificial_domination_gate')
    def test_original_Bs_floor_and_linear_frame_beta_are_kept(self):
        self.control('original_inverse_bound_uses_Bs_floor_squared');self.control('frame_linear_beta_is_not_phase_squared_beta')

    def test_coarse_old_floor_does_not_certify_stronger_phase_inverse_constant(self):self.control('coarse_Bs_floor_does_not_prove_stronger_phase_beta_inverse_bound')
    def test_no_actual_profile_or_stronger_phase_scope_promotion(self):self.control('no_stronger_phase_beta_or_actual_profile_promotion')
    def test_verifier_never_invokes_the_band_selector(self):
        f=fixture();old=selector.select_tail
        def forbidden(*a,**k):raise AssertionError('Band reselection is forbidden')
        selector.select_tail=forbidden
        try:self.assertTrue(verify_existing(f['frame'],f['selector'],f['selected']['ell_tail'])['conditional_normal_floor_implication_passes'])
        finally:selector.select_tail=old

    def test_incomplete_actual_input_sets_and_floats_rejected(self):
        p=frame.inputs();p.pop('Mz')
        with self.assertRaises(ValueError):normal_length_certificate(p)
        p=frame.inputs();p['u']=0.5
        with self.assertRaises(TypeError):normal_length_certificate(p)
        f=fixture();sp=dict(f['selector']);sp.pop('qstar')
        with self.assertRaises(ValueError):verify_existing(f['frame'],sp,16)

    def test_invalid_actual_length_and_slow_errors_rejected(self):
        for S,L,e,e2,ki in ((1,2,0,0,0),(4,0,0,0,0),(4,1,-1,0,0)):
            with self.assertRaises(ValueError):tangential_error(frame.inputs(),S,L,e,e2,ki)

    def test_five_frozen_records_no_added_assumptions_and_source_bindings(self):
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertEqual(len(self.record['dependency_evidence_sha256']),5);self.assertTrue(all(self.record['gates'].values()))

    def test_new_frozen_evidence_reproduces_exactly(self):self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)
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
