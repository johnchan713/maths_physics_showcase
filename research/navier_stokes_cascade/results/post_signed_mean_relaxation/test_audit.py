#!/usr/bin/env python3
"""Current mean maps, surviving slow correlations, pressure and scope controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,generate,validate,compare_records,main
from relax import Harmonic,Poly,R,Z,EPS,KAPPA,avg_product,physical_increment_powers


class MeanRelaxationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.record=generate()

    def review(self,name): return self.record['reviews'][name]

    def test_temporal_inverse_exact_on_both_quadratures(self):
        f=Harmonic(0,R+Z,R*Z)
        self.assertEqual(f.inverse_nt().nt(),f)
        self.assertEqual(f.nt().inverse_nt(),f)

    def test_temporal_inverse_rejects_nonzero_Haar_mean(self):
        with self.assertRaises(ValueError): Harmonic(1,R,Z).inverse_nt()

    def test_exact_product_average_keeps_both_old_new_correlations(self):
        a=Harmonic(R,Z,R*Z);b=Harmonic(Z,R,R+Z)
        self.assertEqual(avg_product(a,b),R*Z+(Z*R+R*Z*(R+Z))/2)
        with self.assertRaises(TypeError): a*b

    def test_graph_average_has_no_added_graph_loss(self):
        f=Harmonic(R*Z,R**2,Z)
        self.assertEqual(f.dr().mean,Z)
        self.assertEqual(f.nt().mean,0)
        self.assertEqual(f.dr().nt(),f.nt().dr())

    def test_fast_actual_graph_curl_is_divergence_free(self):
        self.assertTrue(self.review('fast_inverse_and_nonzero_old_new_correlations')['checks']['graph_curl_divergence'])

    def test_fast_velocity_has_zero_Haar_mean_but_old_slow_mean_does_not(self):
        c=self.review('fast_inverse_and_nonzero_old_new_correlations')['checks']
        self.assertTrue(c['exact_zero_Haar_fast_velocity'])
        self.assertTrue(c['inherited_slow_means_are_not_zero'])
        self.assertTrue(c['fast_update_preserves_existing_slow_mean'])

    def test_fast_source_cancellation_includes_both_desired_components(self):
        self.assertTrue(self.review('fast_inverse_and_nonzero_old_new_correlations')['checks']['temporal_inverse_cancels_both_desired_sources'])

    def test_exact_fast_radial_and_all_three_defect_changes(self):
        c=self.review('fast_inverse_and_nonzero_old_new_correlations')['checks']
        self.assertTrue(all(c[n] for n in ('exact_fast_averaged_radial_source_change',
            'exact_fast_pressure_defect_change','exact_fast_angular_defect_change','exact_fast_axial_defect_change')))

    def test_old_fast_new_fast_correlations_do_not_average_to_zero(self):
        self.assertTrue(self.review('fast_inverse_and_nonzero_old_new_correlations')['checks']['old_fast_new_fast_correlations_survive'])

    def test_stale_prefast_targets_leave_nonzero_defect_changes(self):
        self.assertTrue(self.review('fast_inverse_and_nonzero_old_new_correlations')['checks']['all_three_fast_defect_changes_are_nonzero'])

    def test_slow_radial_remainder_keeps_existing_Haar_mean(self):
        c=self.review('slow_recomputed_defects_with_existing_Haar_mean')['checks']
        self.assertTrue(c['exact_slow_radial_remainder'])
        self.assertTrue(c['dropping_existing_slow_mean_changes_remainder'])

    def test_zero_Haar_parts_drop_only_against_slow_increments(self):
        self.assertTrue(self.review('slow_recomputed_defects_with_existing_Haar_mean')['checks']['postfast_zero_mean_parts_drop_only_from_slow_crosses'])

    def test_slow_curl_and_two_velocity_moments(self):
        c=self.review('slow_recomputed_defects_with_existing_Haar_mean')['checks']
        self.assertTrue(c['slow_graph_curl_divergence'])
        self.assertTrue(c['both_slow_velocity_moments_zero'])

    def test_all_slow_defect_differences_are_recomputed(self):
        c=self.review('slow_recomputed_defects_with_existing_Haar_mean')['checks']
        self.assertTrue(all(c[n] for n in ('exact_pressure_defect_difference','exact_angular_defect_difference','exact_axial_defect_difference')))

    def test_old_slow_cross_terms_have_nonzero_defect_moments(self):
        c=self.review('slow_recomputed_defects_with_existing_Haar_mean')['checks']
        self.assertTrue(c['old_slow_swirl_cross_has_nonzero_radial_moment'])
        self.assertTrue(c['old_slow_flux_crosses_are_nonzero'])
        self.assertTrue(c['nonlinear_recomputed_defects_not_set_to_zero'])

    def test_three_pressures_and_integrated_flux_identities(self):
        r=self.review('all_three_pressures_and_integrated_fluxes')
        self.assertEqual(r['exact_pressure_and_integrated_checks'],15)
        self.assertTrue(r['checks']['all_three_reconstructed_pressure_and_integrated_identities'])

    def test_six_velocity_moment_checks(self):
        r=self.review('all_three_pressures_and_integrated_fluxes')
        self.assertEqual(r['exact_velocity_moment_checks'],6)
        self.assertTrue(r['checks']['all_velocity_moments_preserved'])

    def test_complete_fixed_wave_covariance_and_axial_pressure_retained(self):
        c=self.review('all_three_pressures_and_integrated_fluxes')['checks']
        self.assertTrue(c['all_six_fixed_updated_covariances_retained'])
        self.assertTrue(c['axial_pressure_gradient_not_zero'])
        self.assertTrue(c['moving_bump_derivative_nonzero'])

    def test_pressure_cannot_be_reused_after_slow_update(self):
        self.assertTrue(self.review('all_three_pressures_and_integrated_fluxes')['checks']['stale_fast_pressure_invalid_for_final_state'])

    def test_actual_torus_inverse_not_replaced_by_one_circle_fixture(self):
        self.assertTrue(self.review('fast_inverse_and_nonzero_old_new_correlations')['manufactured_one_circle_graph_not_actual_irrational_two_torus'])
        self.assertTrue(self.review('fast_inverse_and_nonzero_old_new_correlations')['actual_A1_flat_remainder_retained_analytically_and_in_replayed_parent'])

    def test_fast_viscosity_limits_full_mean_gain(self):
        r=self.review('all_fixed_order_mean_gain_and_physical_ledger')
        self.assertEqual(Q(r['fast_remainder_orders']['radial_viscosity']),3-6*KAPPA)
        self.assertEqual(Q(r['rows']['final_full_tangential']),3-6*KAPPA)

    def test_slow_scalar_gain_is_not_twice_target_order(self):
        r=self.review('all_fixed_order_mean_gain_and_physical_ledger');alpha=2-3*KAPPA
        self.assertEqual(Q(r['rows']['slow_updated_defects']),alpha+1)
        self.assertLess(alpha+1,2*alpha)

    def test_all_fast_and_slow_remainder_rows_meet_final_target(self):
        r=self.review('all_fixed_order_mean_gain_and_physical_ledger')
        self.assertGreaterEqual(min(map(Q,r['fast_remainder_orders'].values())),3-6*KAPPA)
        self.assertGreater(min(map(Q,r['slow_tangential_orders'].values())),3-6*KAPPA)

    def test_new_mean_wave_errors_do_not_cancel_old_nonzero_residual(self):
        r=self.review('all_fixed_order_mean_gain_and_physical_ledger')
        self.assertEqual(Q(r['rows']['added_nonzero_wave_source']),2-4*KAPPA)
        self.assertFalse(FLAGS['current_nonzero_angular_residual_cancelled'])

    def test_uniform_physical_velocity_pressure_and_partial_potential_powers(self):
        t=2-4*KAPPA;p=physical_increment_powers(Q(1,200),t,0,0)
        self.assertEqual(p,dict(velocity=Q(-2475001,5000000),pressure=Q(-5000001,5000000),radial_axial_curl_potential=Q(24999,5000000)))
        for a,b in ((1,0),(0,1),(3,2),(6,4)):
            shift=-a*Q(201,400)-b*Q(403,400)
            self.assertEqual(physical_increment_powers(Q(1,200),t,a,b),{n:v+shift for n,v in p.items()})

    def test_physical_parameter_and_derivative_validation(self):
        for h,t,a,b in ((0,2,0,0),(Q(1,100),2,0,0),(Q(1,200),-1,0,0),(Q(1,200),2,-1,0),(Q(1,200),2,True,0)):
            with self.assertRaises(ValueError): physical_increment_powers(h,t,a,b)

    def test_normalized_mean_gain_is_not_absolute_physical_residual_decay(self):
        h=Q(1,200);nu=3-6*KAPPA
        self.assertEqual(h*nu-2*(Q(1,2)+h)-Q(1,2),Q(-14950003,10000000))

    def test_parent_and_625_inherited_hashes_replayed(self):
        self.assertTrue(self.record['gates']['actual_signed_parent_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),625)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_evidence_reproduces(self):
        frozen=json.loads((HERE/'evidence.json').read_text());validate(frozen);compare_records(self.record,frozen)

    def test_false_zero_mean_full_residual_and_blowup_promotions_rejected(self):
        for name,value in FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][name]=True
                with self.assertRaises(ValueError): validate(r)

    def test_current_update_order_and_analytical_dependencies_cannot_change(self):
        for section,name,value in (('protocol','operation_order','reuse old targets'),
            ('flags','relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas',False)):
            r=deepcopy(self.record);r[section][name]=value
            with self.assertRaises(ValueError): validate(r)

    def test_missing_gate_or_obligation_rejected(self):
        for section in ('gates','remaining_obligations'):
            r=deepcopy(self.record)
            if section=='gates': r[section].pop(next(iter(r[section])))
            else: r[section].pop()
            with self.assertRaises(ValueError): validate(r)

    def test_changed_ledger_or_provenance_rejected(self):
        r=deepcopy(self.record);r['reviews']['all_fixed_order_mean_gain_and_physical_ledger']['rows']['slow_updated_defects']='4'
        with self.assertRaises(ValueError): validate(r)
        for section in ('source_hashes','inherited_provenance'):
            r=deepcopy(self.record);r[section][next(iter(r[section]))]=False
            with self.assertRaises(ValueError): validate(r)

    def test_cli_does_not_overwrite_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit): main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__': unittest.main()
