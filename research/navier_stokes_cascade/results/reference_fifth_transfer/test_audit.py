"""Exact higher reference transfer diagnostics and strict scope guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import reference_bounds as r

class ReferenceFifthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.review=cls.record['reviews']
    def test_axis_sixth_coefficient_list(self):self.assertEqual(r.axis_coefficients(6),[1,4,10,20,35,56,84])
    def test_axis_whole_C6_sum(self):self.assertEqual(5*sum(r.axis_coefficients(6)),1050)
    def test_axis_fifth_whole_norm_is_630(self):self.assertEqual(5*sum(r.axis_coefficients(5)),630)
    def test_axis_zero_order_mass(self):self.assertEqual(5*sum(r.axis_coefficients(0)),5)
    def test_invalid_axis_order_rejected(self):
        for k in (-1,True,6.):
            with self.assertRaises(ValueError):r.axis_coefficients(k)
    def test_tame_sum_at_one(self):self.assertLess(r.exponential_majorant(1,6),3)
    def test_tame_sum_at_larger_envelope(self):self.assertLess(r.exponential_majorant(7,6),3*7**6)
    def test_tame_invalid_inputs_rejected(self):
        for K,k in ((0,6),(1.,6),(1,7),(1,True)):
            with self.assertRaises(ValueError):r.exponential_majorant(K,k)
    def test_positive_averaging_contracts(self):self.assertLessEqual(abs(r.convex_average([Q(1,3),Q(2,3)],[-7,5])),7)
    def test_signed_or_wrong_mass_weights_not_accepted(self):
        for w,v in (([-1,2],[1,2]),([1,1],[1,2]),([],[]),([1],[1,2])):
            with self.assertRaises(ValueError):r.convex_average(w,v)
    def test_reference_all_seven_angular_average_identities(self):self.assertTrue(self.review['algebra']['checks']['exact_reference_taper_identity_all_seven_angular_jets'])
    def test_reference_weight_mass_and_sign(self):
        a=self.review['algebra']['checks'];self.assertTrue(a['reference_weights_mass_one']);self.assertTrue(a['reference_weight_nonnegative_by_factorization'])
    def test_eta_dependent_cutoff_not_silently_allowed(self):self.assertTrue(self.review['algebra']['checks']['eta_dependent_taper_commutator_is_not_zero'])
    def test_axis_atoms_points_and_orders_exercised(self):self.assertEqual(self.review['algebra']['axis_diagnostic_coefficients'],112)
    def test_native_sixth_packet_constants_charged(self):
        a=self.review['algebra']['checks'];self.assertTrue(a['native_log_constant_fits_Nnat']);self.assertTrue(a['native_inverse_constant_fits_Nnat']);self.assertTrue(a['native_mixed_log_constant_fits_Nnat'])
    def test_derivative_cost_six_is_charged_in_W(self):self.assertTrue(self.review['algebra']['checks']['W_C5_constant_covers_d_M_eta'])
    def test_complete_source_constants_charged(self):
        a=self.review['algebra']['checks'];self.assertTrue(a['Sq_C5_constant_covers_all_sources']);self.assertTrue(a['Sn_C5_constant_covers_full_pressure'])
    def test_both_integrated_initial_terms_charged(self):
        a=self.review['algebra']['checks'];self.assertTrue(a['p1_C5_initial_and_source_constants_fit']);self.assertTrue(a['ns_C5_initial_and_source_constants_fit'])
    def test_actual_interval_precisions(self):self.assertEqual([v['interval_digits'] for v in self.review['parameters']],[80,110])
    def test_native_C6_packet_comparison(self):self.assertTrue(all(v['checks']['native_packet_below_K'] for v in self.review['parameters']))
    def test_regular_F_and_E_C6_small(self):self.assertTrue(all(v['checks']['reference_F_C6_below_one_over_100'] and v['checks']['reference_E_C6_below_one'] for v in self.review['parameters']))
    def test_complete_reference_pressure_C6_comparison(self):self.assertTrue(all(v['checks']['full_reference_pressure_C6_below_K'] for v in self.review['parameters']))
    def test_actual_core_entrance_five_moment_bounds(self):self.assertTrue(all(v['checks']['five_core_entry_moments_C5_below_K'] for v in self.review['parameters']))
    def test_reference_p1_ns_C5_fit_unchanged_B(self):self.assertTrue(all(v['checks']['reference_p1_C5_below_B_over100'] and v['checks']['reference_ns_C5_below_B_over100'] for v in self.review['parameters']))
    def test_no_enormous_scales_materialized(self):self.assertTrue(all(v['enormous_C_G_H_N_or_actual_profile_materialized'] is False for v in self.review['parameters']))
    def test_endpoints_and_interior_in_source_diagnostic(self):self.assertEqual([v['eta'] for v in self.review['flows']['rows']],['-1','1/3','1'])
    def test_true_forward_moments_satisfy_p1_fifth_ODE(self):self.assertTrue(all(v['p1_angular5_exact'] for v in self.review['flows']['rows']))
    def test_true_forward_moments_satisfy_ns_fifth_ODE(self):self.assertTrue(all(v['ns_angular5_exact'] for v in self.review['flows']['rows']))
    def test_all_mixed_total5_source_coefficients(self):self.assertTrue(all(v['full_source_mixed_total5_exact'] for v in self.review['flows']['rows']))
    def test_dropping_axis_pressure_fifth_detected(self):self.assertTrue(all(v['omitted_axis_fifth_source_nonzero'] for v in self.review['flows']['rows']))
    def test_dropping_radial_pressure_fifth_detected(self):self.assertTrue(all(v['omitted_radial_pressure_fifth_source_nonzero'] for v in self.review['flows']['rows']))
    def test_sixth_trace_and_fifth_fields_nonzero(self):self.assertTrue(all(v['sixth_pressure_trace_nonzero'] and v['fifth_fields_nonzero'] for v in self.review['flows']['rows']))
    def test_manufactured_diagnostics_not_selected_profile(self):self.assertIs(self.review['flows']['manufactured_not_selected_profile'],True)
    def test_protocol_parent_and_manuscript_SHA(self):self.assertEqual(audit.protocol()['parent_commit'],'7a282097f781d13dc410e64d606109aff676f6c2');self.assertEqual(len(audit.protocol()['manuscript_sha256']),64)
    def test_historical_source_bindings(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_record_exact_reproduction(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                q=deepcopy(self.record);q['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(q)
    def test_numeric_zero_not_boolean_false(self):
        q=deepcopy(self.record);q['flags']['independent_foundation_review_or_formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(q)
    def test_domain_axis_memory_taper_or_scope_not_erased(self):
        for k in ('domain','axis','reference','tame','sources','transfer','entry','scope'):
            q=deepcopy(self.record);q['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(q)
    def test_next_actual_continuation_obligation_retained(self):
        q=deepcopy(self.record);q['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(q)
    def test_fresh_evidence_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(path)]);audit.validate(json.loads(path.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(path)])

if __name__=='__main__':unittest.main()
