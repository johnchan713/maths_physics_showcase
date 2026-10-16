#!/usr/bin/env python3
"""Current accumulated means, complete complex covariance, pressure and scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,generate,validate,compare_records,main
from relaxation import KAPPA,Harmonic,R,Z,T,fixture,avg_product,physical_powers


class ComplexMeanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.record=generate()

    def review(self,name): return self.record['reviews'][name]
    def check(self,review,name): self.assertTrue(self.review(review)['checks'][name])

    def test_temporal_inverse_on_both_current_quadratures(self):
        f=Harmonic(0,R+Z,R*T)
        self.assertEqual(f.inverse_nt().nt(),f);self.assertEqual(f.nt().inverse_nt(),f)

    def test_temporal_inverse_rejects_nonzero_Haar_source(self):
        with self.assertRaises(ValueError): Harmonic(1,R,Z).inverse_nt()

    def test_accumulated_Haar_mean_includes_previous_slow_update(self):
        self.check('current_fast_inverse_and_accumulated_Haar_state','accumulated_Haar_mean_contains_previous_slow_update')

    def test_full_current_complex_tensor_differs_from_old_real_signed_tensor(self):
        self.check('complete_fixed_complex_tensor_pressure_and_velocity_moments','actual_fixture_covariance_updated_beyond_first_real_signed_state')

    def test_fast_actual_increment_has_zero_Haar_mean(self):
        self.check('current_fast_inverse_and_accumulated_Haar_state','zero_Haar_actual_fast_increment')

    def test_fast_graph_curl_is_divergence_free(self):
        self.check('current_fast_inverse_and_accumulated_Haar_state','complete_graph_curl_divergence')

    def test_current_accumulated_Haar_mean_is_preserved_by_fast_map(self):
        self.check('current_fast_inverse_and_accumulated_Haar_state','preserves_accumulated_old_slow_Haar_mean')

    def test_fast_inverse_cancels_both_current_sources(self):
        self.check('current_fast_inverse_and_accumulated_Haar_state','temporal_inverse_cancels_both_current_sources')

    def test_exact_fast_radial_source_change(self):
        self.check('current_fast_inverse_and_accumulated_Haar_state','exact_fast_radial_source_change')

    def test_all_fast_defect_changes_recomputed(self):
        for n in ('exact_fast_pressure_defect_change','exact_fast_angular_defect_change','exact_fast_axial_defect_change'):
            self.check('current_fast_inverse_and_accumulated_Haar_state',n)

    def test_old_fast_correlations_and_stale_target_errors_nonzero(self):
        for n in ('old_new_fast_correlation_is_nonzero','stale_prefast_targets_leave_three_nonzero_changes'):
            self.check('current_fast_inverse_and_accumulated_Haar_state',n)

    def test_exact_slow_radial_remainder_uses_accumulated_means(self):
        self.check('current_slow_defects_and_stale_target_obstruction','exact_slow_radial_remainder_with_accumulated_means')

    def test_slow_zero_mean_simplification_only_after_Haar_average(self):
        self.check('current_slow_defects_and_stale_target_obstruction','only_zero_Haar_parts_drop_against_slow_increment')
        self.check('current_slow_defects_and_stale_target_obstruction','deleting_old_Haar_mean_changes_remainder')

    def test_slow_curl_and_both_velocity_moments(self):
        for n in ('exact_slow_graph_curl_divergence','two_slow_velocity_moments_zero'):
            self.check('current_slow_defects_and_stale_target_obstruction',n)

    def test_current_five_map_inputs_not_numerically_replaced_by_fixture(self):
        r=self.review('current_slow_defects_and_stale_target_obstruction')
        self.assertTrue(r['actual_fixed_five_map_cancels_current_postfast_linear_rows_analytically'])
        self.assertTrue(r['manufactured_delta_not_numerical_actual_five_map_solution'])

    def test_all_slow_defect_differences_recomputed(self):
        for n in ('exact_slow_pressure_defect_difference','exact_slow_angular_defect_difference','exact_slow_axial_defect_difference'):
            self.check('current_slow_defects_and_stale_target_obstruction',n)

    def test_surviving_old_slow_cross_terms_and_nonlinear_defects_nonzero(self):
        for n in ('old_swirl_cross_has_nonzero_radial_moment','old_axial_and_swirl_cross_fluxes_nonzero','nonlinear_recomputed_defects_are_not_zero'):
            self.check('current_slow_defects_and_stale_target_obstruction',n)

    def test_three_full_pressure_and_integrated_flux_reconstructions(self):
        r=self.review('complete_fixed_complex_tensor_pressure_and_velocity_moments')
        self.assertEqual(r['exact_pressure_flux_checks'],15)
        self.assertTrue(r['checks']['pressure_and_integrated_fluxes_at_all_three_current_states'])

    def test_six_velocity_moment_checks(self):
        r=self.review('complete_fixed_complex_tensor_pressure_and_velocity_moments')
        self.assertEqual(r['exact_velocity_moment_checks'],6)
        self.assertTrue(r['checks']['all_six_velocity_moments_zero'])

    def test_axial_pressure_and_whole_moving_bump_not_zero(self):
        for n in ('axial_mean_pressure_gradient_nonzero','moving_bump_derivative_nonzero'):
            self.check('complete_fixed_complex_tensor_pressure_and_velocity_moments',n)

    def test_stale_pressure_invalid_after_both_maps(self):
        self.check('complete_fixed_complex_tensor_pressure_and_velocity_moments','pressure_changes_after_both_mean_operations')

    def test_current_full_tangential_order_comes_from_averaged_covariance(self):
        r=self.review('current_mean_gain_physical_and_uncancelled_wave_ledger')
        nu=2-6*KAPPA
        self.assertEqual(Q(r['rows']['final_full_tangential']),nu)
        self.assertGreater(min(map(Q,r['fast_remainder_orders'].values())),nu)
        self.assertGreater(min(map(Q,r['slow_tangential_orders'].values())),nu)

    def test_defects_gain_alpha_plus_one_not_twice_alpha(self):
        r=self.review('current_mean_gain_physical_and_uncancelled_wave_ledger')['rows'];alpha=Q(r['alpha'])
        self.assertEqual(Q(r['slow_updated_defects']),alpha+1)
        self.assertLess(alpha+1,2*alpha)

    def test_precise_stale_target_shortfall_is_two_kappa(self):
        r=self.review('current_mean_gain_physical_and_uncancelled_wave_ledger')['rows']
        self.assertEqual(Q(r['slow_updated_defects'])-Q(r['fast_defect_change']),2*KAPPA)

    def test_old_nonzero_residual_is_unchanged_and_limiting(self):
        r=self.review('current_mean_gain_physical_and_uncancelled_wave_ledger')['rows']
        self.assertEqual(Q(r['final_nonzero_regular']),Q(3,2)-6*KAPPA)
        self.assertLess(Q(r['final_nonzero_regular']),Q(r['final_full_tangential']))

    def test_added_mean_wave_sources_and_current_fixed_primary_remainder(self):
        r=self.review('current_mean_gain_physical_and_uncancelled_wave_ledger')['rows']
        self.assertEqual(Q(r['added_nonzero_wave_source']),Q(3,2)-4*KAPPA)
        self.assertEqual(Q(r['current_fixed_primary_remainder']),1-3*KAPPA)

    def test_all_fixed_physical_derivative_powers(self):
        h=Q(1,200);H=Q(3,2)-4*KAPPA;p=physical_powers(h,H,0,0)
        self.assertEqual(p,dict(velocity=Q(-2487501,5000000),pressure=Q(-5012501,5000000),radial_axial_curl_potential=Q(12499,5000000)))
        for a,b in ((1,0),(0,1),(3,2),(6,4)):
            shift=-a*Q(201,400)-b*Q(403,400)
            self.assertEqual(physical_powers(h,H,a,b),{n:v+shift for n,v in p.items()})

    def test_physical_parameter_and_derivative_validation(self):
        for h,t,a,b in ((0,1,0,0),(Q(1,100),1,0,0),(Q(1,200),-1,0,0),(Q(1,200),1,-1,0),(Q(1,200),1,True,0)):
            with self.assertRaises(ValueError): physical_powers(h,t,a,b)

    def test_normalized_gain_not_absolute_physical_residual_decay(self):
        h=Q(1,200);D=Q(3,2)-6*KAPPA
        self.assertEqual(-2*(Q(1,2)+h)-Q(1,2)+h*D,Q(-15025003,10000000))

    def test_manufactured_graph_not_actual_global_two_torus_inverse(self):
        self.assertTrue(self.review('current_fast_inverse_and_accumulated_Haar_state')['one_circle_manufactured_inverse_not_actual_global_two_torus_solve'])

    def test_parent_and_671_inherited_hashes(self):
        self.assertTrue(self.record['gates']['current_generic_signed_parent_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),671)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_evidence_reproduces(self):
        frozen=json.loads((HERE/'evidence.json').read_text());validate(frozen);compare_records(self.record,frozen)

    def test_false_scope_promotions_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError): validate(r)

    def test_changed_operation_order_or_analytical_dependencies_rejected(self):
        for section,n,v in (('protocol','operation_order','reuse stale targets'),
            ('flags','relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas',False)):
            r=deepcopy(self.record);r[section][n]=v
            with self.assertRaises(ValueError): validate(r)

    def test_changed_ledger_missing_gate_obligation_or_provenance_rejected(self):
        r=deepcopy(self.record);r['reviews']['current_mean_gain_physical_and_uncancelled_wave_ledger']['rows']['final_nonzero_regular']='2'
        with self.assertRaises(ValueError): validate(r)
        for section in ('gates','remaining_obligations','source_hashes','inherited_provenance'):
            r=deepcopy(self.record)
            if section=='remaining_obligations': r[section].pop()
            else: r[section].pop(next(iter(r[section])))
            with self.assertRaises(ValueError): validate(r)

    def test_cli_does_not_overwrite_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit): main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__': unittest.main()
