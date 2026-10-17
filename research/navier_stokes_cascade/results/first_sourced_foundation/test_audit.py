#!/usr/bin/env python3
"""Complete first sourced-step identities and strict scope controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from sourced_bounds import KAPPA,orders,physical_powers,inverse_polynomial_losses
from sourced_checks import forced_checks,tensor_checks,pressure_checks
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main


class FirstSourcedFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_live_parent_and_unchanged_operations(self):
        self.assertEqual(protocol()['parent_commit'],'9dbaa3cf768f3009f4f17d9b5c5b47e4df613a6b')
        self.assertIn('operation order',protocol()['unchanged_parameters'])

    def test_signed_and_mean_prefix_preserves_actual_B(self):
        o=orders();self.assertEqual(o['B'],1-KAPPA);self.assertTrue(all(x>=o['B'] for x in o['prefix'].values()))

    def test_first_sourced_complete_regular_groups_and_gain(self):
        o=orders();self.assertEqual(min(o['generated'].values()),F(3,2)-2*KAPPA)
        self.assertEqual(o['source_gain'],F(1,2)-KAPPA)

    def test_actual_complete_covariance_change_C_includes_square(self):
        o=orders();self.assertEqual(o['covariance_change'],o['C']);self.assertGreaterEqual(2*o['B'],o['C'])

    def test_unaveraged_and_Haar_pressure_orders_are_distinct(self):
        o=orders();self.assertEqual(o['averaged_mean_pressure_change']-o['full_mean_pressure_change'],KAPPA)

    def test_current_defects_and_full_means_not_frozen_at_pre_wave_orders(self):
        o=orders();self.assertLess(o['recomputed_defects'],o['pre_source_defects']);self.assertLess(o['new_full_tangential_mean'],o['pre_source_full_mean'])

    def test_next_old_wave_error_is_B_not_residual_D(self):
        o=orders();self.assertEqual(o['next_generic_old_error'],o['B']);self.assertLess(o['old_error'],o['D'])

    def test_all_fixed_Duhamel_loss_degrees_remain_finite(self):
        g,b=inverse_polynomial_losses([2*j+1 for j in range(15)])
        self.assertEqual(g,[2*j+1 for j in range(15)]);self.assertEqual(b,[2*j+2 for j in range(15)])

    def test_bad_loss_degrees_and_kappa_rejected(self):
        for v in ([],[-1],[True]):
            with self.assertRaises(ValueError):inverse_polynomial_losses(v)
        for k in (0,F(1,10),-1):
            with self.assertRaises(ValueError):orders(k)

    def test_all_physical_derivative_shifts_match_fixed_phase_costs(self):
        h=F(1,200);p0=physical_powers(h,0,0)
        for a in range(5):
            for b in range(4):
                p=physical_powers(h,a,b);loss=a*(F(1,2)+h/2)+b*(1+3*h/2)
                self.assertTrue(all(p[k]==v-loss for k,v in p0.items()))

    def test_normalized_gain_is_not_physical_absolute_residual_decay(self):self.assertLess(physical_powers(F(1,200),0,0)['regular_residual'],0)

    def test_invalid_physical_h_or_jet_count_rejected(self):
        for h,a,b in ((0,0,0),(F(1,100),0,0),(F(1,200),True,0),(F(1,200),0,-1)):
            with self.assertRaises(ValueError):physical_powers(h,a,b)

    def test_all_600_exact_complex_forced_pressure_and_constraint_components(self):
        r=forced_checks();self.assertEqual(r['exact_mixed_constraint_pressure_ODE_components'],600);self.assertTrue(r['checks']['complex_moving_constraint_and_forced_pressure_exact'])

    def test_both_source_normal_and_moving_normal_pressure_terms_nonzero(self):
        c=forced_checks()['checks'];self.assertTrue(c['source_normal_pressure_term_nonzero']);self.assertTrue(c['moving_normal_pressure_term_nonzero'])

    def test_negative_harmonics_preserve_complex_reality(self):self.assertTrue(forced_checks()['checks']['conjugate_negative_harmonics_preserve_reality'])
    def test_zero_entry_does_not_give_temporal_zero_exit(self):self.assertTrue(forced_checks()['checks']['zero_entry_does_not_force_zero_exit'])

    def test_independent_Fourier_zero_modes_match_all_36_components(self):
        r=tensor_checks();self.assertEqual(r['exact_zero_mode_components'],36);self.assertTrue(r['checks']['independent_complex_Fourier_zero_mode_matches_complete_tensor'])

    def test_complete_old_and_sourced_curls_remain_divergence_free(self):self.assertTrue(tensor_checks()['checks']['complete_old_and_new_curls_divergence_free'])
    def test_all_nine_complete_tensor_increment_components(self):
        r=tensor_checks();self.assertEqual(r['exact_tensor_increment_components'],9);self.assertTrue(r['checks']['full_cross_plus_all_new_squares_identity'])

    def test_sourced_sine_principal_breaks_initial_real_shortcut(self):self.assertTrue(tensor_checks()['checks']['old_sine_new_principal_sine_cross_nonzero'])
    def test_second_harmonic_square_not_lost_to_orthogonality(self):self.assertTrue(tensor_checks()['checks']['second_harmonic_square_not_dropped'])
    def test_unequal_harmonics_not_merged_before_angular_average(self):self.assertTrue(tensor_checks()['checks']['unequal_angular_harmonics_must_not_be_merged'])
    def test_generic_covariance_change_need_not_be_PSD(self):self.assertTrue(tensor_checks()['checks']['generic_covariance_increment_need_not_be_positive_semidefinite'])

    def test_full_before_after_pressure_recomputed_with_nonzero_means(self):
        r=pressure_checks();self.assertEqual(r['exact_before_after_pressure_and_moment_controls'],14);self.assertTrue(all(r['checks'].values()))

    def test_whole_moving_pressure_bump_derivative_retained(self):
        c=pressure_checks()['checks'];self.assertTrue(c['whole_moving_pressure_product']);self.assertTrue(c['moving_bump_derivative_omission_detected'])

    def test_same_full_background_pressure_not_leading_truncation(self):
        c=pressure_checks()['checks'];self.assertTrue(c['same_full_background_pressure_subtraction']);self.assertTrue(c['truncating_background_pressure_would_leave_tail_gradient'])

    def test_exact_velocity_moments_not_changed_by_wave_operation(self):self.assertTrue(pressure_checks()['checks']['exact_old_velocity_moments_unchanged'])

    def test_finite_polynomial_regularity_not_C_infinity_certificate(self):
        self.assertIn('C4 only',tensor_checks()['fixture_potential_zero_extension'])
        self.assertTrue(forced_checks()['fixture_not_actual_giant_pulse_or_endpoint_source_certificate'])

    def test_no_all_stage_residual_PSD_or_velocity_promotion(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)

    def test_all_907_inherited_and_18_direct_source_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),907);self.assertEqual(len(self.record['source_hashes']),18)

    def test_frozen_record_and_JSON_roundtrip_match(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r);validate(json.loads(json.dumps(self.record)))

    def test_cached_exact_reviews_not_mutable_evidence_aliases(self):
        r=generate();r['reviews']['complex_forced_constraint_and_pressure']['checks']['complex_moving_constraint_and_forced_pressure_exact']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(forced_checks()['checks']['complex_moving_constraint_and_forced_pressure_exact'])

    def test_same_count_provenance_key_substitution_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])

    def test_protocol_gates_sources_and_remaining_obligations_not_erased(self):
        for key in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[key],dict):r[key].pop(next(iter(r[key])))
            else:r[key].pop()
            with self.assertRaises(ValueError):validate(r)

    def test_cli_refuses_overwriting_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')


if __name__=='__main__':unittest.main()
