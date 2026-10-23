#!/usr/bin/env python3
"""Meaningful exact normalization mutations and stage/pressure checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import sympy as s
import audit
from controls import review
from assembly import principal_defect,rates


class LocalizedStressTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        values=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(values,fragment);self.assertTrue(all(values),fragment)

    def test_independent_cutoff_integrals_and_physical_Haar_measure(self):self.control('two_cutoffs_angular_half_and_Haar_J_retained')
    def test_positive_original_roots_in_manufactured_patch(self):self.control('positive_roots_on_manufactured_patch')
    def test_H_selected_flux_agrees_with_full_tensor_integrals(self):self.control('cutoff_columns_match_direct_tensor_integrals')
    def test_two_different_bands_have_same_physical_stress(self):self.control('physical_scale_normalization_exact')
    def test_squared_partition_and_both_sign_convention(self):self.control('squared_partition_counted_once_per_box')
    def test_same_localized_principal_stress_is_exact(self):self.control('full_assembled_selected_stress_exact')
    def test_all_six_entries_remain(self):self.control('all_six_covariance_entries_retained')
    def test_three_term_defect_identity_with_simultaneous_errors(self):
        self.control('unconditional_three_defect_identity');self.control('simultaneous_three_nonzero_defects_retained')
    def test_entire_identity_survives_spatial_and_scale_derivatives(self):self.control('entire_partition_identity_survives_')
    def test_nonconstant_partition_derivatives_are_not_dropped_individually(self):self.control('partition_derivatives_nonzero_before_exact_sum')
    def test_nonunit_pulse_cutoff_does_not_create_an_extra_remainder(self):self.control('unmatched_cutoff_not_needed_for_matched_H')
    def test_recutted_pulse_requires_new_column_or_explicit_defect(self):self.control('detect_stale_roots_after_recutted_pulse')
    def test_missing_pulse_cutoff_changes_covariance(self):
        self.control('detect_missing_temporal_cutoff');self.control('detect_missing_transverse_cutoff')
    def test_missing_angular_half_or_Haar_J_is_detected(self):
        self.control('detect_missing_angular_half');self.control('detect_missing_Haar_J')
    def test_wrong_square_velocity_scale_is_detected(self):self.control('detect_wrong_physical_velocity_square_scale')
    def test_duplicate_signs_or_missing_boxes_are_detected(self):
        self.control('detect_signs_counted_twice_as_boxes');self.control('detect_incomplete_active_cover')
    def test_negative_controls_still_obey_unconditional_defect_identity(self):self.control('defect_identity_for_')
    def test_full_postwave_tangential_order_D_for_entire_kappa_range(self):self.control('full_tangential_postwave_minimum_is_D')
    def test_stage_is_old_means_with_new_wave(self):self.control('postwave_uses_initial_old_means_not_fresh_final')
    def test_complete_rows_retain_leading_gap_until_substitution(self):
        self.control('leading_defect_retained_in_both_complete_rows');self.control('combines_Q1_and_entire_Haar_rows')
    def test_pressure_compensator_and_original_high_jet_A1_remain(self):
        self.control('full_pressure_source_and_compensator_retained');self.control('exact_original_A1_high_jet_recipe_retained')
    def test_tails_and_actual_input_gate_not_promoted(self):
        self.control('additive_tails_excluded');self.control('no_actual_input_certificate');self.control('no_full_radial_tail_or_actual_input_promotion')
    def test_final_Haar_bounds_preserved(self):self.control('final_Haar_and_radial_Haar_preserved')
    def test_missing_labels_and_incomplete_box_data_rejected(self):
        with self.assertRaises(ValueError):principal_defect(s.Matrix([1,2]),[])
        with self.assertRaises(ValueError):principal_defect(s.Matrix([1,2]),[dict(eta=1)])
    def test_invalid_kappa_rejected(self):
        for k in (0,Q(1,4),-1,0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
    def test_frozen_evidence_and_six_dependencies(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),6)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_remaining_obligations_and_dependency_mutations_rejected(self):
        for key,value in self.record['flags'].items():
            bad=deepcopy(self.record);bad['flags'][key]=not value
            with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['dependency_evidence_sha256'][next(iter(bad['dependency_evidence_sha256']))]='wrong'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_existing_frozen_record_not_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
