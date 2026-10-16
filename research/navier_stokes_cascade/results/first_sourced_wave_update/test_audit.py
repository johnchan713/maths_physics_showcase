#!/usr/bin/env python3
"""First sourced inverse, flat-tail obstruction, complex pressure and scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,generate,validate,compare_records,main
from sources import KAPPA,physical_powers


class FirstSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.record=generate()

    def review(self,name): return self.record['reviews'][name]
    def check(self,review,name): self.assertTrue(self.review(review)['checks'][name])

    def test_full_cutoff_commutator_includes_viscosity_and_pressure(self):
        self.check('full_viscous_potential_pressure_cutoff','full_linear_potential_pressure_cutoff_commutator')

    def test_cutoff_curl_product_identity(self):
        self.check('full_viscous_potential_pressure_cutoff','exact_curl_product_decomposition')

    def test_cut_potential_velocity_is_divergence_free(self):
        self.check('full_viscous_potential_pressure_cutoff','cut_potential_is_divergence_free')

    def test_all_four_omitted_commutator_contributions_detected(self):
        r=self.review('full_viscous_potential_pressure_cutoff')
        self.assertEqual(set(r['omitted_term_gaps']),{'curl','pressure','viscosity','transport'})
        self.assertTrue(all(v>1e-6 for v in r['omitted_term_gaps'].values()))

    def test_forced_projected_and_constraint_equations_exact(self):
        self.check('constrained_sourced_pressure_and_endpoint','forced_constraint_pressure_and_projected_ODE_exact')
        self.assertEqual(self.review('constrained_sourced_pressure_and_endpoint')['exact_rational_checks'],64)

    def test_pressure_needs_moving_normal(self):
        self.check('constrained_sourced_pressure_and_endpoint','moving_normal_term_is_needed')

    def test_pressure_needs_source_normal_component(self):
        self.check('constrained_sourced_pressure_and_endpoint','source_normal_component_is_not_dropped')

    def test_pressure_needs_imaginary_factor_and_both_harmonic_signs(self):
        for n in ('missing_pressure_i_is_detected','both_positive_and_negative_harmonics_checked'):
            self.check('constrained_sourced_pressure_and_endpoint',n)

    def test_zero_entry_does_not_imply_zero_exit(self):
        self.check('constrained_sourced_pressure_and_endpoint','zero_entry_does_not_imply_zero_exit')

    def test_early_step_derivative_integral_and_maximum(self):
        for n in ('normalized_early_step_derivative_integral_one','exact_derivative_sup_majorant'):
            self.check('exact_forward_inverse_flat_tail_obstruction',n)

    def test_finite_step_regularities_not_conflated(self):
        r=self.review('exact_forward_inverse_flat_tail_obstruction')
        self.assertIn('C2',r['diagnostic_step_regularity']);self.assertIn('C-infinity',r['diagnostic_step_regularity'])
        self.check('exact_forward_inverse_flat_tail_obstruction','polynomial_step_has_zero_first_two_endpoint_jets')

    def test_gaussian_small_source_response_stays_nonflat(self):
        r=self.review('exact_forward_inverse_flat_tail_obstruction')
        self.assertTrue(all(x['inverse_core_value']==-1 for x in r['rows']))
        self.assertLess(r['rows'][-1]['log_upper_source_norm'],-96)
        self.check('exact_forward_inverse_flat_tail_obstruction','Gaussian_source_norm_bound_decreases')

    def test_source_zero_endpoint_collars_still_does_not_suffice(self):
        r=self.review('exact_forward_inverse_flat_tail_obstruction')
        self.assertTrue(all(x['source_zero_on_both_endpoint_collars'] for x in r['rows']))
        self.check('exact_forward_inverse_flat_tail_obstruction','physical_flatness_is_not_preserved_by_forward_inverse')

    def test_tail_is_retained_not_erased_or_inverted(self):
        r=self.review('exact_forward_inverse_flat_tail_obstruction')
        self.assertTrue(r['tail_kept_as_residual_not_erased_or_forward_inverted'])
        self.assertEqual(r['exact_actual_homogeneous_witness'],'f=chi_prime*t_h; solution=(-chi*t_h,-chi*pi_h)')

    def test_general_two_quadrature_covariance_keeps_crosses_and_square(self):
        self.check('general_complex_covariance_and_axial_pressure','general_two_quadrature_full_tensor_increment')
        self.assertEqual(self.review('general_complex_covariance_and_axial_pressure')['exact_tensor_checks'],18)

    def test_initial_real_shortcut_is_not_reused(self):
        self.check('general_complex_covariance_and_axial_pressure','initial_real_shortcut_not_reused')

    def test_all_six_covariance_changes_nonzero(self):
        self.check('general_complex_covariance_and_axial_pressure','all_six_covariance_changes_nonzero')

    def test_compact_full_axial_pressure_moment(self):
        for n in ('compact_pressure_reconstructed','axial_pressure_moment','changed_axial_pressure_and_moving_bump_nonzero'):
            self.check('general_complex_covariance_and_axial_pressure',n)

    def test_both_integrated_fluxes_and_all_three_defects(self):
        for n in ('angular_integrated_flux','axial_integrated_flux','all_three_defects_must_be_recomputed'):
            self.check('general_complex_covariance_and_axial_pressure',n)

    def test_current_source_minimum_and_low_order_tail(self):
        r=self.review('regular_source_output_and_physical_ledger')
        self.assertEqual(min(map(Q,r['current_source_terms'].values())),1-3*KAPPA)
        self.assertEqual(Q(r['rows']['homogeneous_principal_tail_order']),Q(1,2))

    def test_complete_new_products_and_viscous_linear_minimum(self):
        r=self.review('regular_source_output_and_physical_ledger')
        self.assertEqual(set(r['new_source_terms']),{'full_linear','old_new_wave','new_wave_square','actual_mean_cross'})
        self.assertEqual(min(map(Q,r['new_source_terms'].values())),Q(3,2)-6*KAPPA)

    def test_covariance_means_and_defects_regenerated(self):
        r=self.review('regular_source_output_and_physical_ledger')['rows']
        for n in ('complete_covariance_change','averaged_mean_pressure_change','new_Haar_tangential_mean','new_defects'):
            self.assertEqual(Q(r[n]),Q(3,2)-3*KAPPA)
        self.assertEqual(Q(r['new_full_tangential_mean']),Q(3,2)-4*KAPPA)

    def test_next_fixed_primary_remainder_is_lower_than_old_curl(self):
        r=self.review('regular_source_output_and_physical_ledger')['rows']
        self.assertEqual(Q(r['old_fixed_primary_remainder_after_solve']),1-3*KAPPA)
        self.assertLess(Q(r['old_fixed_primary_remainder_after_solve']),1-KAPPA)

    def test_one_finite_joined_order_gain(self):
        r=self.review('regular_source_output_and_physical_ledger')['rows']
        self.assertEqual(Q(r['finite_joined_order']),Q(3,2)-6*KAPPA)
        self.assertEqual(Q(r['finite_order_gain']),Q(1,2)-3*KAPPA)

    def test_finite_harmonics_support_iteration_not_uniformity(self):
        r=self.review('regular_source_output_and_physical_ledger')
        self.assertEqual(r['current_harmonic_multipliers'],[-2,-1,1,2])
        self.assertEqual(r['maximum_updated_residual_multiplier'],4)
        self.assertTrue(r['finite_stage_constants_not_uniform_in_iteration_or_a_fixed_scale_contraction'])

    def test_all_fixed_physical_derivative_powers(self):
        h=Q(1,200);B=1-3*KAPPA;p=physical_powers(B,h,0,0)
        self.assertEqual(p,dict(velocity=Q(-10000003,20000000),pressure=Q(-20050003,20000000),potential=Q(49997,20000000)))
        for a,b in ((1,0),(0,1),(3,2),(6,4)):
            shift=-a*Q(201,400)-b*Q(403,400)
            self.assertEqual(physical_powers(B,h,a,b),{n:v+shift for n,v in p.items()})

    def test_physical_parameter_validation(self):
        for t,h,a,b in ((-1,Q(1,200),0,0),(1,0,0,0),(1,Q(1,100),0,0),(1,Q(1,200),-1,0),(1,Q(1,200),True,0)):
            with self.assertRaises(ValueError): physical_powers(t,h,a,b)

    def test_normalized_gain_does_not_imply_absolute_physical_decay(self):
        h=Q(1,200);D=Q(3,2)-6*KAPPA
        self.assertEqual(-2*(Q(1,2)+h)-Q(1,2)+h*D,Q(-15025003,10000000))

    def test_diagnostics_not_actual_global_fields(self):
        self.assertTrue(self.review('full_viscous_potential_pressure_cutoff')['manufactured_physical_fields_not_actual_pulses'])
        self.assertTrue(self.review('general_complex_covariance_and_axial_pressure')['arbitrary_complex_quadrature_fixture_not_an_actual_sourced_pulse'])

    def test_parent_and_640_inherited_hashes(self):
        self.assertTrue(self.record['gates']['current_post_signed_parent_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),640)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_evidence_reproduces(self):
        frozen=json.loads((HERE/'evidence.json').read_text());validate(frozen);compare_records(self.record,frozen)

    def test_false_scope_promotions_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError): validate(r)

    def test_tail_split_and_analytical_dependencies_cannot_change(self):
        for section,n,v in (('protocol','exact_source_split','invert everything'),
            ('flags','relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas',False)):
            r=deepcopy(self.record);r[section][n]=v
            with self.assertRaises(ValueError): validate(r)

    def test_changed_ledger_missing_gate_and_provenance_rejected(self):
        r=deepcopy(self.record);r['reviews']['regular_source_output_and_physical_ledger']['rows']['D']='2'
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
