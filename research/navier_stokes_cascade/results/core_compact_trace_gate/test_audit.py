"""Trace defects and negative controls, distinct from canonical validation."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from bounds import boundary_recipe,seam_orders,A1_trace_request,certify_actual_traces,transported_modes
from controls import review

class TraceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_weighted_characteristic_drift_for_all_three_e(self):self.control('weighted_drift_pullback_retains_full_phase')
    def test_projection_seam_sign(self):self.control('whole_classical_derivative_has_negative_internal_jump');self.control('original_chi_R_and_sign')
    def test_cumulative_jump_and_same_cutoff_in_primitive(self):self.control('cumulative_jump_and_cutoff_repair')
    def test_total_zero_jump_is_not_enough_for_primitive(self):self.control('zero_total_jump_does_not_remove_cumulative_primitive_defect')
    def test_pressure_local_primitive_needs_matching(self):self.control('local_pressure_not_equal_to_primitive_without_trace_matching')
    def test_value_gluing_recovers_smooth_identity(self):self.control('smooth_value_gluing_restores_exact_primitive')
    def test_moving_cutoff_seam_pressure_derivative(self):self.control('moving_cutoff_derivative_of_seam_pressure_is_retained')
    def test_zero_Haar_moments_do_not_certify_traces(self):self.control('two_Haar_moments_do_not_detect')
    def test_repeated_IBP_includes_all_seams_and_exterior(self):self.control('all_seam_and_exterior_traces_reconstruct_exact_integral')
    def test_bulk_only_IBP_is_detected(self):self.control('bulk_only_mutation_detects_boundary_terms')
    def test_finite_regular_compact_leak_matches_independent_convolution(self):self.control('finite_regular_compact_fixture_has_exact_M_minus_n_leak')
    def test_fixed_mode_leak_is_nonzero(self):self.control('nonzero_leak_has_exact_magnitude')
    def test_lower_jets_of_comparison_really_match(self):self.control('all_lower_derivative_traces_match')
    def test_moving_integral_has_interface_speed_term(self):self.control('moving_integral_keeps_seam_speed_times_value_jump')
    def test_dropping_moving_seam_changes_derivative(self):self.control('dropping_moving_seam_term_changes_classical_derivative')
    def test_mixed_trace_chain_rule(self):self.control('mixed_trace_chain_rule_keeps_next_radial_trace')
    def test_continuity_of_flux_does_not_imply_continuity_of_slow_derivative(self):self.control('continuous_Phi_alone_does_not_make_Dz_Phi_value_trace_zero')
    def test_exact_matching_and_last_trace_budget(self):self.control('all_exact_matching_recovers_bulk_power');self.control('last_trace_is_allowed_at_bulk_order')
    def test_unmatched_trace_blocks_separate_norm_argument(self):self.control('smallest_unmatched_trace_blocks_strict_target')
    def test_extra_quantified_trace_smallness_can_suffice(self):self.control('quantified_extra_epsilon_smallness_can_replace_exact_matching')
    def test_original_physical_requests_not_lowered(self):self.control('physical_p30000001');self.control('Phi_through30000000')
    def test_direct_normalized_route_needs_no_new_IBP(self):self.control('normalized_direct_A1_has_no_added')
    def test_giant_recipe_deferred(self):self.control('giant_trace_request_is_deferred')
    def test_actual_input_certification_refused(self):
        self.control('no_actual_canonical_trace_or_complete_decay_claim')
        with self.assertRaisesRegex(NotImplementedError,'canonical whole closed-cover'):certify_actual_traces()
    def test_incomplete_trace_mapping_rejected(self):
        with self.assertRaisesRegex(ValueError,'every trace'):seam_orders(2,3,{0:None})
    def test_equality_is_not_strict_margin(self):
        r=seam_orders(2,2,{0:Q(0),1:None},kappa=Q(1,10),target=Q(21,10))
        self.assertFalse(r['strict_target_pass'])
    def test_invalid_integer_orders_and_floats_rejected(self):
        for n in (-1,True,Q(1,2),1.0):
            with self.assertRaises(ValueError):boundary_recipe(n)
        with self.assertRaises(ValueError):boundary_recipe(0)
        with self.assertRaises(TypeError):seam_orders(2.0,1,{0:None})
    def test_boundary_recipe_is_constant_size(self):
        self.assertEqual(len(boundary_recipe(10)),len(boundary_recipe(100000000)))
    def test_physical_trace_uses_original_cover_transport(self):
        self.control('original_J_and_nonunimodular_cover_retained');self.control('matching_physical_trace_uses_transposed_relative_cover')
        self.control('equal_raw_mode_lists_are_not_equal_physical_traces');self.control('relative_cover_composition_uses_same_original_matrix')
    def test_Haar_preservation_does_not_prove_cover_trace_matching(self):
        self.control('zero_Haar_is_preserved_but_does_not_make_traces_match');self.control('constant_mode_Haar_is_unchanged')
    def test_invalid_or_unbounded_relative_cover_rejected(self):
        with self.assertRaises(ValueError):transported_modes({(1,0):Q(1)},4)
        with self.assertRaises(ValueError):transported_modes({(Q(1,2),0):Q(1)},1)
    def test_all_controls_and_six_records_pass(self):
        self.assertTrue(all(self.r['checks'].values()))
        self.assertTrue(all(self.record['gates'].values()))
        self.assertEqual(len(self.record['dependency_evidence_sha256']),6)
    def test_record_round_trip(self):audit.validate(json.loads(json.dumps(self.record)))
    def test_record_rejects_scope_pressure_or_trace_mutations(self):
        for key in ('canonical_closed_cover_traces_or_all_scale_norms_verified','full_current_pressure_memory_and_axial_viscosity_retained'):
            r=deepcopy(self.record);r['flags'][key]=not r['flags'][key]
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['review']['details']['budget']['certificate']['physical_illustration']['integrations']=30000000
        with self.assertRaises(ValueError):audit.validate(r)
    def test_frozen_record(self):
        p=Path(__file__).with_name('evidence.json')
        if p.exists():audit.validate(json.loads(p.read_text()))
    def test_audit_refuses_overwrite_and_reproduces_fresh_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'fresh.json';audit.main(['--output',str(p)])
            audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
