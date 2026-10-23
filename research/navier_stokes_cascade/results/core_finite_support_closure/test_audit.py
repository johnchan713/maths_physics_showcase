"""Finite-support closure factors and mutation checks, not actual norm samples."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review
from bounds import flat_radius_recipe,signed_support_threshold,certify_actual_inputs,support_recipe

class SupportClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_all_mixed_signed_jets_keep_half_weight(self):self.control('every_signed_mixed_jet_retains_half_weight')
    def test_root_derivatives_cannot_be_dropped(self):self.control('root_denominator_derivatives_cannot_be_dropped')
    def test_original_root_product_only_is_extended(self):self.control('original_root_product_not_reciprocal_alone_is_extended')
    def test_sharp_target_exponential_weight_threshold(self):self.control('sharp_support_threshold_has_explicit_noncanonical_comparison')
    def test_smooth_target_alone_does_not_supply_signed_support(self):self.control('smooth_compact_half_weight_target_can_give_nonvanishing_signed_coefficient')
    def test_invalid_target_weight_rejected(self):
        for theta in (0,-1):
            with self.assertRaises(ValueError):signed_support_threshold(theta)
        with self.assertRaises(TypeError):signed_support_threshold(0.5)
    def test_flat_edge_jets_have_only_finite_poles(self):self.control('flat_edge_derivative_is_exponential_times_finite_poles');self.control('no_negative_exponential_weight_created')
    def test_moving_edge_mixed_chain_rule(self):self.control('moving_chart_full_RZT_chain_matches_direct_flat_factor')
    def test_arbitrary_fixed_poles_absorbed_with_explicit_strict_power(self):self.control('arbitrary_finite_poles_leave_strict_vanishing_power')
    def test_zero_entrance_full_jets_give_zero_response_jets(self):self.control('zero_full_entrance_jets_force_zero_response_jets')
    def test_Gaussian_smallness_is_not_entrance_flatness(self):self.control('Gaussian_smallness_cannot_replace_entrance_flatness')
    def test_raw_outgoing_response_is_not_assumed_compact(self):self.control('homogeneous_outgoing_response_is_not_automatically_compact')
    def test_whole_path_face_support_is_preserved(self):self.control('whole_path_parameter_face_support_is_preserved')
    def test_entrance_vanishing_does_not_supply_path_face_support(self):self.control('entrance_vanishing_alone_does_not_supply_whole_path_face_support')
    def test_complete_cylindrical_cut_potential_curl(self):self.control('complete_cut_potential_keeps_cutoff_curl_and_connection')
    def test_velocity_only_cutoff_rejected(self):self.control('velocity_only_cutoff_mutation_changes_field')
    def test_full_cut_pressure_gradient(self):self.control('full_cut_pressure_slow_derivative_keeps_cutoff_term')
    def test_source_normal_pressure_remains_when_velocity_zero(self):self.control('zero_projected_forcing_velocity_does_not_remove_full_source_normal_pressure')
    def test_compact_inverse_can_have_nonzero_input_moment(self):self.control('nonzero_source_moment_is_allowed_by_compact_Te')
    def test_potential_hull_fills_gap(self):self.control('actual_potential_hull_includes_gap_between_source_and_cutoff')
    def test_same_cutoff_restores_exterior_potential(self):self.control('same_cutoff_restores_zero_exterior_potential')
    def test_global_physical_representative_not_inverse_branch(self):self.control('global_integer_pullback_has_one_physical_phase');self.control('periodic_cover_uses_no_independent_inverse_branch')
    def test_conditional_traces_follow_without_mass_sampling(self):self.control('all_mixed_trace_orders_follow_conditionally_not_from_samples')
    def test_old_jet_graph_and_physical_request_retained(self):self.control('same_old_physical_raw_request_and_entire_DAG_retained')
    def test_actual_input_certification_refused(self):
        self.control('actual_weighted_input_or_same_band_gaps_not_certified')
        with self.assertRaisesRegex(NotImplementedError,'Actual weighted mixed'):certify_actual_inputs()
    def test_large_vanishing_constant_not_materialized(self):
        a=flat_radius_recipe(60000006,30000000)
        self.assertFalse(a['materialized']);self.assertGreater(a['resulting_power'],30000000)
    def test_invalid_fixed_orders_rejected(self):
        for n in (-1,True,1.0):
            with self.assertRaises(ValueError):flat_radius_recipe(n,0)
            with self.assertRaises(ValueError):support_recipe(n)
    def test_all_controls_and_six_records(self):
        self.assertTrue(all(self.r['checks'].values()));self.assertTrue(all(self.record['gates'].values()))
        self.assertEqual(len(self.record['dependency_evidence_sha256']),6)
    def test_round_trip_and_frozen_evidence(self):
        audit.validate(json.loads(json.dumps(self.record)))
        p=Path(__file__).with_name('evidence.json')
        if p.exists():audit.validate(json.loads(p.read_text()))
    def test_scope_and_pressure_mutations_rejected(self):
        for key in ('original_weighted_canonical_input_or_same_band_floor_gaps_certified','full_cut_potential_curl_and_source_normal_pressure_retained','negative_physical_time_zero_extension_assumed'):
            r=deepcopy(self.record);r['flags'][key]=not r['flags'][key]
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['protocol']['means_pressure']='discard pressure'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--output',str(p)])
            audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
