"""Axis regularity, full-pressure parity, Cauchy constants and conditional scope."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import axis_bounds as b

class AxisFirstTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_every_A0_entry_has_reflection_covariance(self):self.assertTrue(self.r['symbolic']['checks']['all_A0_entries_have_reflection_covariance'])
    def test_every_A1_entry_has_reflection_covariance(self):self.assertTrue(self.r['symbolic']['checks']['all_A1_entries_have_reflection_covariance'])
    def test_pressure_and_both_axial_sources_reflect(self):self.assertTrue(self.r['symbolic']['checks']['full_pressure_and_axial_source_has_reflection_covariance'])
    def test_diagonal_commutes_with_reflection(self):self.assertTrue(self.r['symbolic']['checks']['diagonal_commutes_with_reflection'])
    def test_variable_kernel_block_annihilation(self):self.assertTrue(self.r['symbolic']['checks']['variable_diagonal_kernel_derivative_blocks_annihilate'])
    def test_derivative_of_right_block_annihilation(self):self.assertTrue(self.r['symbolic']['checks']['derivative_hitting_right_block_also_annihilates'])
    def test_distinct_radial_kernels_annihilation(self):self.assertTrue(self.r['symbolic']['checks']['different_radii_kernels_do_not_spoil_annihilation'])
    def test_complex_ray_Volterra_monomials(self):self.assertTrue(self.r['symbolic']['checks']['complex_ray_regular_Volterra_monomial_identities'])
    def test_odd_source_parity_negative_control(self):self.assertTrue(self.r['symbolic']['checks']['odd_swirl_source_breaks_needed_parity'])
    def test_72_exact_scalar_coefficient_recurrences(self):self.assertEqual(self.r['recurrence']['exact_scalar_recurrences'],72);self.assertTrue(self.r['recurrence']['checks']['all_six_component_recurrences'])
    def test_all_12_orders_have_component_parities(self):self.assertEqual(self.r['recurrence']['terms'],12);self.assertTrue(self.r['recurrence']['checks']['all_component_parities'])
    def test_zero_axis_traces(self):self.assertTrue(self.r['recurrence']['checks']['zero_axis_data'])
    def test_regular_average_component_identity(self):self.assertTrue(self.r['recurrence']['checks']['exact_regular_average_component'])
    def test_only_first_four_components_descend_to_X(self):self.assertTrue(self.r['recurrence']['checks']['first_four_components_descend_to_holomorphic_X']);self.assertTrue(self.r['recurrence']['checks']['last_two_are_odd_not_declared_holomorphic_X'])
    def test_pressure_omission_changes_axial_coefficient(self):self.assertTrue(self.r['recurrence']['checks']['removing_pressure_source_changes_axial_coefficient'])
    def test_parity_break_creates_odd_profile(self):self.assertTrue(self.r['recurrence']['checks']['parity_breaking_source_creates_odd_regular_profile'])
    def test_actual_kinematic_first_axis_slopes_in_manufactured_system(self):self.assertEqual(self.r['recurrence']['first_axis_X_slopes'],['e/8 + 1/8','1/2 - e/4','-e**2/2 - 1/2'])
    def test_complex_core_majorant_includes_radial_order_two(self):self.assertTrue(self.r['envelopes']['checks']['core_majorant_bounds_all_needed_radial_orders'])
    def test_four_delta_within_core_joint_margin(self):self.assertTrue(self.r['envelopes']['checks']['actual_four_delta_strip_inside_core_joint_majorant'])
    def test_complex_core_Cauchy_disk_margin(self):self.assertTrue(self.r['envelopes']['checks']['Cauchy_disks_stay_strictly_inside_complex_core'])
    def test_core_and_positive_collar_overlap_at_all_scales(self):
        for L in (F(1,1000),F(1),F(4),F(1000)):
            j=b.joined_jets(L,F(1,4),2,3);self.assertTrue(j['overlap']);self.assertLessEqual(j['xi_positive_lower']**2,j['Xa'])
    def test_core_radial_Cauchy_factors(self):self.assertEqual(b.core_envelope(2,1,3)['radial_factors'],(1,6,72,1296))
    def test_join_keeps_two_distinct_M_constants(self):
        j=b.joined_jets(3,F(1,4),2,3);self.assertEqual(j['radial_Ceta2_bounds'][0],(73,73,0));self.assertEqual(j['radial_Ceta2_bounds'][3],(94608,F(68401,8),F(2555,8)))
    def test_common_bound_dominates_every_jet_component(self):
        j=b.joined_jets(3,F(1,4),2,3)
        self.assertTrue(all(j['common_mixed_bound'][i]>=row[i] for row in j['radial_Ceta2_bounds'] for i in range(3)))
        self.assertEqual(j['common_raw_mixed_bound'],tuple(2*x for x in j['common_mixed_bound']))
    def test_second_angular_Cauchy_cost_not_first_only(self):self.assertEqual(b.joined_jets(3,F(1,4),2,3)['Ceta2_factor'],73)
    def test_complex_path_length_bound(self):self.assertTrue(self.r['envelopes']['checks']['complex_path_length_bound_valid'])
    def test_full_source_retains_shifted_viscosities(self):self.assertTrue(self.r['envelopes']['checks']['full_pressure_source_envelope_dominates_both_shifted_viscosities'])
    def test_exponential_remains_symbolic(self):self.assertTrue(self.r['envelopes']['checks']['no_selected_exponentials_materialized'])
    def test_Cauchy_monomial_controls(self):self.assertTrue(self.r['envelopes']['checks']['analytic_monomials_obey_core_Cauchy_factors'])
    def test_third_cutoff_jet_remains_a_distinct_input(self):self.assertTrue(self.r['envelopes']['checks']['third_cutoff_jet_not_discharged_by_axis_join'])
    def test_cut_product_identity_cutoff(self):
        rows=b.joined_jets(3,F(1,4),2,3)['radial_Ceta2_bounds'];self.assertEqual(b.cut_product(rows,(1,0,0,0)),rows)
    def test_three_original_matrix_structure_samples(self):self.assertEqual(len(self.r['original_structure']['rows']),3);self.assertTrue(self.r['original_structure']['checks']['original_system_structure_matches'])
    def test_invalid_core_inputs(self):
        for S,rho,L in ((1,1,1),(2,0,1),(2,2,1),(2,1,0)):
            with self.assertRaises(ValueError):b.core_envelope(S,rho,L)
    def test_invalid_join_inputs(self):
        for L,delta,Q,Qxi in ((0,F(1,4),2,3),(1,0,2,3),(1,2,2,3),(1,F(1,4),0,3),(1,F(1,4),2,-1)):
            with self.assertRaises(ValueError):b.joined_jets(L,delta,Q,Qxi)
    def test_exact_inputs_reject_float_bool(self):
        for bad in (2.,True):
            with self.assertRaises(TypeError):b.core_envelope(bad,1,1)
    def test_cut_product_bad_shapes_and_signs(self):
        with self.assertRaises(ValueError):b.cut_product(((1,2,3),)*3,(1,0,0,0))
        with self.assertRaises(ValueError):b.cut_product(((1,2,3),)*4,(1,0,0,-1))
    def test_invalid_recurrence_order(self):
        for n in (0,-1,True):
            with self.assertRaises(ValueError):b.recurrence(n)
    def test_parent_protocol_and_frozen_record(self):self.assertEqual(audit.protocol()['parent_commit'],'d93abfc5528b2b6b4826144973ba5aeed574eed9');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_sources_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_unsupported_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_not_false(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_core_activation_and_pressure_scope_not_erased(self):
        for k in ('core','cutoff','parity'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_second_H23_obligation_not_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(1)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_frozen_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
