"""Seventh input, actual sixth cover, same nonconstant root and native R5 scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import native6 as n
class NativeSixthTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
 def test_native_seventh_generating_factor(self):self.assertTrue(self.r['algebra']['checks']['seventh_coefficient_geometric_factor_below_ten'])
 def test_native_sixth_n_factor18(self):self.assertTrue(self.r['algebra']['checks']['sixth_differentiated_n_factor_below_eighteen'])
 def test_native_log_inverse_and_mixed_constants(self):self.assertTrue(self.r['algebra']['checks']['native_log_inverse_C7_constants_fit_N7']);self.assertTrue(self.r['algebra']['checks']['native_mixed_C6_constant_fits_N7'])
 def test_axis_seventh_positive_mixture1650(self):self.assertEqual(self.r['algebra']['axis_C7_coefficient_sum'],1650)
 def test_positive_taper_all_eight_jets(self):self.assertTrue(self.r['algebra']['checks']['positive_taper_all_eight_angular_jets_exact'])
 def test_positive_mass_one_weights(self):self.assertTrue(self.r['algebra']['checks']['positive_taper_mass_one']);self.assertTrue(self.r['algebra']['checks']['positive_taper_weight_nonnegative'])
 def test_eta_dependent_taper_not_substituted(self):self.assertTrue(self.r['algebra']['checks']['eta_dependent_taper_would_not_commute'])
 def test_seventh_finite_exp_and_T7_constant(self):self.assertTrue(self.r['algebra']['checks']['seventh_finite_exponential_sum_below_three']);self.assertTrue(self.r['algebra']['checks']['T7_constant_covers_all_seven_exp_jets'])
 def test_reference_sixth_derivative_cost7(self):self.assertTrue(self.r['algebra']['checks']['sixth_angular_derivative_cost_is_seven'])
 def test_reference_sources_keep_both_pressure_terms(self):self.assertTrue(self.r['algebra']['checks']['reference_Sq6_constant_covers_all_terms']);self.assertTrue(self.r['algebra']['checks']['reference_Sn6_constant_covers_both_pressure_terms'])
 def test_reference_variation_constants(self):self.assertTrue(self.r['algebra']['checks']['reference_p1_C6_variation_constant_fits']);self.assertTrue(self.r['algebra']['checks']['reference_ns_C6_variation_constant_fits'])
 def test_actual_shape_log_and_exp_sixth(self):self.assertTrue(self.r['algebra']['checks']['sixth_shape_log_f_still_below_six']);self.assertTrue(self.r['algebra']['checks']['sixth_shape_exp_sum_still_below1000'])
 def test_ideal_f_sixth_majorant28(self):self.assertTrue(self.r['algebra']['checks']['sixth_f_majorant_sum28'])
 def test_normalizer_C6_still5_25(self):self.assertEqual(self.r['algebra']['normalizer_C6_costs'],['5','25'])
 def test_annular_energy_constant8125(self):self.assertTrue(self.r['algebra']['checks']['annular_energy_constant_still8125'])
 def test_same_eta_dependent_sixth_root_identity(self):self.assertTrue(self.r['root']['checks']['sixth_full_eta_dependent_root_recurrence_exact'])
 def test_freezing_earlier_Q_detected(self):self.assertTrue(self.r['root']['checks']['freezing_earlier_Q_eta_detected_at_sixth'])
 def test_Q4_c1_c1_retained(self):self.assertTrue(self.r['root']['checks']['Q4_c1_c1_term_present_at_sixth'])
 def test_no_need_for_Q5_Q6(self):self.assertTrue(self.r['root']['checks']['no_unneeded_fifth_sixth_Q_derivatives'])
 def test_seventh_reference_inputs_sixth_source_identity(self):self.assertTrue(self.r['flows']['checks']['seventh_input_sixth_reference_source_identity'])
 def test_axis_seventh_nonzero(self):self.assertTrue(self.r['flows']['checks']['seventh_axis_trace_retained'])
 def test_deleting_axis_or_Cp_detected(self):self.assertTrue(self.r['flows']['checks']['deleting_axis_sixth_source_detected']);self.assertTrue(self.r['flows']['checks']['deleting_Cp_sixth_source_detected'])
 def test_all_three_actual_stage_shears(self):self.assertTrue(self.r['flows']['checks']['actual_all_three_prescription_stage_shears_total6_exact']);self.assertEqual(len(self.r['flows']['actual_stage_rows']),9)
 def test_all_five_actual_densities(self):self.assertTrue(self.r['flows']['checks']['actual_all_five_density_flows_total6_exact'])
 def test_full_actual_pressure_flow(self):self.assertTrue(self.r['flows']['checks']['actual_full_pressure_flow_total6_exact'])
 def test_sixth_incoming_moments_not_zeroed(self):self.assertTrue(self.r['flows']['checks']['six_order_incoming_moment_traces_not_zeroed'])
 def test_actual_moments_not_reference(self):self.assertTrue(self.r['flows']['checks']['actual_moments_not_copied_from_reference'])
 def test_exact_fixture_endpoints_and_interior(self):self.assertEqual([r['eta'] for r in self.r['flows']['reference_rows']],['-1','1/3','1'])
 def test_native_E_U_regional_SUM(self):self.assertTrue(self.r['envelope']['checks']['sum_of_all_native_E_C6_regions_fits_A2']);self.assertTrue(self.r['envelope']['checks']['sum_of_all_native_U_C6_regions_fits_A'])
 def test_native_one_radial_fifth_SUM(self):self.assertTrue(self.r['envelope']['checks']['native_first_radial_region_SUM_fits_A5']);self.assertEqual(self.r['envelope']['native_first_radial_C5_bound_power'],5)
 def test_full_moment_sixth_budget_A8(self):self.assertTrue(self.r['envelope']['checks']['all_five_moments_plus_regular_Cp_and_axis_fit_A8']);self.assertEqual(self.r['envelope']['full_moment_C6_bound_power'],8)
 def test_inverse_E_fifth_A17(self):self.assertTrue(self.r['envelope']['checks']['inverse_E_C5_fits_A17'])
 def test_shear_fifth_A23(self):self.assertTrue(self.r['envelope']['checks']['native_shear_C5_fits_A23'])
 def test_W_Qs_Ns_keep_full_rows(self):
  for k in ('six_normalized_moment_derivative_cost_W_fits_A10','complete_Ns_with_axial_pressure_fits_A12','complete_Qs_fits_A29'):self.assertTrue(self.r['envelope']['checks'][k])
 def test_complete_both_pressure_C5_A31(self):self.assertTrue(self.r['envelope']['checks']['complete_p1_C5_fits_A31']);self.assertTrue(self.r['envelope']['checks']['complete_p2_C5_fits_A31'])
 def test_native_R5_A32(self):self.assertTrue(self.r['envelope']['checks']['actual_native_R5_common_sum_fits_A32']);self.assertEqual(self.r['envelope']['R5_bound_power'],32)
 def test_auxiliary_envelope_worst_P(self):self.assertTrue(all(n.envelope(2**16,2**16)['checks'].values()))
 def test_auxiliary_invalid_parameters(self):
  for a,p in ((True,16),(2**16,1),(2**16,2**17),(2**16,16.),(1,16)):
   with self.assertRaises(ValueError):n.envelope(a,p)
 def test_interval_precisions80_110(self):self.assertEqual([r['interval_digits'] for r in self.r['parameters']],[80,110]);self.assertTrue(all(all(r['checks'].values()) for r in self.r['parameters']))
 def test_actual_N7_below_unchanged_Kold(self):self.assertTrue(all(r['checks']['native_seventh_packet_below_same_Kold'] for r in self.r['parameters']))
 def test_actual_J6_below_unchanged_Bold(self):self.assertTrue(all(r['checks']['actual_J6_below_same_B_over100'] for r in self.r['parameters']))
 def test_same_C6_root_radius_and_contraction(self):self.assertTrue(all(r['checks']['same_C6_root_contraction'] and r['checks']['same_C6_root_coefficients_below_5p72e14'] for r in self.r['parameters']))
 def test_fixed_first_radial_bump_cost_charged(self):self.assertTrue(all(r['checks']['fixed_width_first_radial_cost_fits_C20'] for r in self.r['parameters']))
 def test_no_actual_cutoff_FINAL_or_huge_parameters_evaluated(self):self.assertTrue(all(r['actual_C_G_H_N_Z5_or_FINAL_profile_materialized'] is False for r in self.r['parameters']))
 def test_interval_invalid_precision(self):
  for d in (30,80.,True):
   with self.assertRaises(ValueError):n.parameters(d)
 def test_parent_and_protocol(self):self.assertEqual(audit.PARENT,'1ebd014108a145266238882a7c3b447946088270');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
 def test_native_R5_independent_J(self):self.assertIs(self.record['flags']['native_R5_estimate_independent_of_variance_cutoff_and_joint_target_J'],True)
 def test_no_complete_J6_or_FINAL_promotion(self):self.assertIs(self.record['flags']['complete_native_joint_J6_or_FINAL_modulated_sixth_field_bound_derived'],False)
 def test_fixed136_historical_records(self):h=self.record['inherited_provenance'];self.assertEqual(len(h['frozen_records']),136);self.assertEqual(h['fixed_manifest_sha256'],audit.MANIFEST_SHA)
 def test_all_historical_bindings(self):self.assertIs(self.record['inherited_provenance']['all_bindings_verified'],True);self.assertGreater(self.record['inherited_provenance']['binding_count'],2436)
 def test_exact_frozen_evidence(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
 def test_cutoff_FINAL_global_scope_promotions_rejected(self):
  for k,v in audit.FLAGS.items():
   if not v:
    r=deepcopy(self.record);r['flags'][k]=True
    with self.assertRaises(ValueError):audit.validate(r)
 def test_zero_cannot_replace_false(self):
  r=deepcopy(self.record);r['flags']['actual_Z4_Z5_or_joint_target_J_or_reduced_fifth_cutoff_gate_verified']=0
  with self.assertRaises(ValueError):audit.validate(r)
 def test_all_protocol_bound(self):
  for k in ('premises','native','reference','actual','root','cover','pressure','gate','scope'):
   r=deepcopy(self.record);r['protocol'][k]='erased'
   with self.assertRaises(ValueError):audit.validate(r)
 def test_actual_cutoff_still_next_obligation(self):
  r=deepcopy(self.record);r['remaining_obligations'].pop(0)
  with self.assertRaises(ValueError):audit.validate(r)
 def test_R5_cannot_be_promoted_to_FINAL_acceptance(self):
  r=deepcopy(self.record);r['reviews']['envelope']['actual_variance_cutoff_or_FINAL_gate_verified']=True
  with self.assertRaises(ValueError):audit.validate(r)
 def test_fresh_replay_and_overwrite_refusal(self):
  with tempfile.TemporaryDirectory() as d:
   out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)]);audit.validate(json.loads(out.read_text()))
   with self.assertRaises(SystemExit):audit.main(['--output',str(out)])
if __name__=='__main__':unittest.main()
