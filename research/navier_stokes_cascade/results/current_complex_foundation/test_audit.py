#!/usr/bin/env python3
"""Current generic complex correction and refreshed mean regression gates."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from complex_join_bounds import KAPPA,orders,generic_gain,physical_powers
from complex_join_checks import wave_review,mean_review,pressure_review
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main


class CurrentComplexFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_live_parent_and_unchanged_operations(self):
        self.assertEqual(protocol()['parent_commit'],'bd707aadc121c83563490495c412ab4c7d15e026')
        self.assertIn('operation order',protocol()['unchanged_parameters'])

    def test_current_inputs_not_old_conservative_orders(self):
        o=orders();self.assertEqual((o['B'],o['C'],o['D']),(1-KAPPA,F(3,2)-KAPPA,F(3,2)-2*KAPPA))

    def test_all_three_generic_covariance_error_groups_tie(self):
        o=orders();self.assertEqual(set(o['covariance'].values()),{2-2*KAPPA})

    def test_generic_strict_gain_derived_not_real_phase_shortcut(self):
        self.assertEqual(orders()['generic_gain'],F(1,2)-KAPPA)
        self.assertEqual(generic_gain(1,1),0);self.assertEqual(generic_gain(2,F(1,2)),0)

    def test_invalid_generic_gain_and_kappa_inputs_rejected(self):
        for a,r,k in ((0,1,KAPPA),(2,0,KAPPA),(2,1,0),(2,1,F(1,2))):
            with self.assertRaises(ValueError):generic_gain(a,r,k)
        for k in (0,-1,F(1,10)):
            with self.assertRaises(ValueError):orders(k)

    def test_full_tensor_change_only_C_not_nu(self):self.assertLess(orders()['C'],orders()['nu'])

    def test_complete_signed_operator_and_interactions_preserve_D(self):
        o=orders();self.assertEqual(min(o['signed_nonzero'].values()),o['D'])

    def test_fast_viscosity_above_averaged_error_by_half_minus_two_kappa(self):
        o=orders();self.assertEqual(min(o['fast'].values())-o['nu'],F(1,2)-2*KAPPA)

    def test_fast_defects_change_but_remain_in_C(self):self.assertGreater(orders()['fast_defect_change'],orders()['C'])

    def test_slow_scalar_gain_uses_accumulated_Haar_order_one(self):
        o=orders();self.assertEqual(o['slow_defects'],o['C']+1);self.assertLess(o['slow_defects'],2*o['C'])

    def test_stale_prefast_targets_miss_new_defect_bound_by_two_kappa(self):self.assertEqual(orders()['stale_target_gap'],2*KAPPA)

    def test_final_FULL_tangential_rows_reach_two_minus_two_kappa(self):self.assertEqual(orders()['final_full_tangential'],2-2*KAPPA)

    def test_fast_mean_wave_source_TIES_D_not_above(self):
        o=orders();self.assertEqual(o['mean_wave_fast'],o['D']);self.assertGreater(o['mean_wave_slow'],o['D'])

    def test_regular_nonzero_residual_and_old_primary_error_not_promoted(self):
        o=orders();self.assertEqual(o['final_nonzero'],o['D']);self.assertEqual(o['old_wave_error'],o['B']);self.assertLess(o['B'],o['D'])

    def test_final_seed_meets_mean_half_order_and_stronger_defects(self):
        o=orders();self.assertEqual(o['nu'],o['D']+F(1,2));self.assertGreater(o['slow_defects'],o['nu'])

    def test_every_fixed_physical_jet_has_expected_phase_loss(self):
        h=F(1,200);o=orders()
        for t,wave in ((o['C'],True),(o['H'],False),(o['C'],False)):
            p0=physical_powers(h,t,0,0,wave)
            for a in range(5):
                for b in range(4):
                    p=physical_powers(h,t,a,b,wave);loss=a*(F(1,2)+h/2)+b*(1+3*h/2)
                    self.assertTrue(all(p[n]==v-loss for n,v in p0.items()))

    def test_invalid_physical_domain_or_jet_rejected(self):
        for h,t,a,b in ((0,1,0,0),(F(1,100),1,0,0),(F(1,200),-1,0,0),(F(1,200),1,True,0),(F(1,200),1,0,-1)):
            with self.assertRaises(ValueError):physical_powers(h,t,a,b)

    def test_exact_independent_complex_Fourier_and_divergence_components(self):
        r=wave_review();self.assertEqual(r['exact_Fourier_components'],36);self.assertEqual(r['exact_divergence_components'],6)
        self.assertTrue(r['checks']['independent_multiharmonic_Fourier_zero_mode']);self.assertTrue(r['checks']['all_complete_curls_divergence_free'])

    def test_full_generic_tensor_and_all_six_changes(self):
        c=wave_review()['checks'];self.assertTrue(c['generic_four_group_full_tensor_identity']);self.assertTrue(c['all_six_tensor_changes_nonzero'])

    def test_all_three_generic_error_groups_and_old_complex_error_nonzero(self):
        c=wave_review()['checks'];self.assertTrue(c['each_generic_error_group_has_nonzero_radial_tangential_entries']);self.assertTrue(c['deleting_old_complex_error_changes_increment'])

    def test_old_second_harmonic_retained_in_current_covariance(self):self.assertTrue(wave_review()['checks']['old_second_harmonic_not_deleted_from_current_covariance'])

    def test_fast_zero_Haar_and_exact_graph_divergence(self):
        c=mean_review()['checks'];self.assertTrue(c['fast_all_three_actual_velocities_zero_Haar']);self.assertTrue(c['fast_and_slow_complete_graph_curl_divergence'])

    def test_accumulated_old_slow_mean_not_replaced_by_initial(self):
        c=mean_review()['checks'];self.assertTrue(c['fast_preserves_accumulated_nonzero_old_Haar_means']);self.assertTrue(c['current_old_Haar_contains_previous_slow_increment'])

    def test_fast_three_defect_changes_and_stale_target_negative_control(self):
        c=mean_review()['checks'];self.assertTrue(c['exact_fast_all_three_defect_changes']);self.assertTrue(c['stale_prefast_targets_leave_three_nonzero_changes'])

    def test_slow_full_remainder_and_three_defect_differences(self):
        c=mean_review()['checks']
        for n in ('exact_slow_remainder_with_all_old_means','exact_slow_pressure_defect_difference','exact_slow_angular_defect_difference','exact_slow_axial_defect_difference'):
            self.assertTrue(c[n])

    def test_old_mean_cross_and_nonlinear_slow_remainders_not_zero(self):
        c=mean_review()['checks'];self.assertTrue(c['old_slow_swirl_cross_cannot_be_deleted']);self.assertTrue(c['nonlinear_slow_defect_remainders_not_zero'])

    def test_four_pressure_states_twenty_flux_controls_eight_moments(self):
        r=pressure_review();self.assertEqual(r['exact_pressure_flux_checks'],20);self.assertEqual(r['exact_velocity_moment_checks'],8);self.assertTrue(all(r['checks'].values()))

    def test_moving_pressure_product_axial_gradient_and_each_pressure_refresh(self):
        c=pressure_review()['checks']
        for n in ('moving_bump_derivative_not_dropped','every_mean_pressure_has_nonzero_axial_gradient','pressure_changes_at_each_operation'):
            self.assertTrue(c[n])

    def test_full_background_tail_subtraction_not_leading_only(self):self.assertTrue(pressure_review()['checks']['leading_only_subtraction_leaves_both_tail_gradients'])

    def test_fixtures_not_actual_C_infinity_fields_or_inverse_solutions(self):
        self.assertIn('C4 only',wave_review()['fixture_potential_zero_extension'])
        self.assertTrue(mean_review()['manufactured_increments_NOT_actual_temporal_inverse_or_five_map_solution'])

    def test_all_925_inherited_and_18_direct_source_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),925);self.assertEqual(len(self.record['source_hashes']),18)

    def test_frozen_record_and_JSON_roundtrip_match(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r);validate(json.loads(json.dumps(self.record)))

    def test_no_all_stage_stale_target_residual_or_blowup_promotion(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)

    def test_cached_reviews_not_mutable_aliases(self):
        r=generate();r['reviews']['current_multiharmonic_signed_tensor']['checks']['generic_four_group_full_tensor_identity']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(wave_review()['checks']['generic_four_group_full_tensor_identity'])

    def test_same_count_provenance_substitution_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])

    def test_protocol_gates_sources_remaining_not_erased(self):
        for key in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[key],dict):r[key].pop(next(iter(r[key])))
            else:r[key].pop()
            with self.assertRaises(ValueError):validate(r)

    def test_cli_refuses_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')


if __name__=='__main__':unittest.main()
