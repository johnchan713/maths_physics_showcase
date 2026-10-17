#!/usr/bin/env python3
"""Current compact-force comparison and independent infinity flux controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from compact_bounds import old,young_decay,separated_radius,kernel,flux_rows,pressure_review,radius_review,limit_review
from compact_checks import matrix,local_review,exterior_review
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main

class ComparisonFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_same_current_force_parent_is_exact(self):self.assertEqual(protocol()['parent_commit'],'44fc3767b90e07f11df62b4726d2d48a3d8e85ec')
    def test_same_force_datum_and_pressure_retained(self):self.assertIn('zero datum',protocol()['current_input']);self.assertTrue(FLAGS['common_force_divergence_cancels_exactly_once_in_difference'])
    def test_competitor_has_no_added_derivative_or_pressure_assumptions(self):self.assertIn('no global pressure derivative or dissipation assumption',protocol()['competitor'])
    def test_recovery_uses_momentum_not_Poisson_alone(self):self.assertIn('time-averaged momentum',protocol()['recovery']);self.assertFalse(FLAGS['pressure_Poisson_equation_alone_used_to_remove_harmonic_gradient'])
    def test_negative_Sobolev_orders_avoid_critical_endpoint(self):self.assertIn('H^-2',protocol()['recovery']);self.assertIn('H^-3',protocol()['recovery']);self.assertGreater(2,F(3,2))
    def test_mixed_stress_keeps_both_orders(self):self.assertIn('w tensor u+u tensor w',protocol()['pressure'])
    def test_exact_stress_decomposition_has_nine_components(self):self.assertEqual(local_review()['exact_stress_components'],9)
    def test_full_local_energy_identity_at_three_viscosities(self):self.assertEqual(local_review()['exact_weighted_local_energy_identities'],3)
    def test_every_energy_term_in_fixture_is_nonzero(self):self.assertEqual(local_review()['nonzero_terms'],['dissipation','force','laplacian','pressure','strain','time','transport'])
    def test_difference_source_divergence_not_dropped_before_same_force(self):self.assertTrue(local_review()['checks']['full_forced_difference_pressure_keeps_minus_div_H'])
    def test_local_diagnostics_not_whole_space_solutions(self):self.assertTrue(local_review()['manufactured_local_polynomials_and_difference_residual_not_whole_space_solutions'])
    def test_support_separation_includes_boundary(self):self.assertEqual(separated_radius(2,1),1);self.assertGreaterEqual(separated_radius(7,2),F(7,2))
    def test_bad_radius_or_support_rejected(self):
        for args in ((F(1,2),0),(1,1),(2,-1)):
            with self.assertRaises(ValueError):separated_radius(*args)
    def test_Hessian_radial_sign_and_trace(self):self.assertEqual(matrix((F(1),F(0),F(0))),[[2,0,0],[0,-1,0],[0,0,-1]])
    def test_nonunit_kernel_direction_rejected(self):
        with self.assertRaises(ValueError):matrix((F(1),F(1),F(0)))
    def test_four_rational_unit_direction_controls(self):self.assertEqual(exterior_review()['exact_kernel_direction_controls'],4)
    def test_sixty_exact_compact_exterior_kernel_controls(self):self.assertEqual(exterior_review()['exact_separated_exterior_controls'],60);self.assertTrue(all(abs(F(r['value_without_inverse_four_pi']))<=F(r['bound_without_inverse_four_pi']) for r in exterior_review()['rows']))
    def test_point_stress_controls_not_actual_L1_fields(self):self.assertTrue(exterior_review()['exterior_point_tensor_kernel_controls_not_actual_L1_stresses_or_flows'])
    def test_compact_mixed_pressure_flux_is_five_halves(self):self.assertTrue(pressure_review()['checks']['cross_flux_decay_is_five_halves']);self.assertIn('R^-5/2',pressure_review()['mixed_flux'])
    def test_no_false_L1_pressure_endpoint(self):self.assertTrue(pressure_review()['checks']['no_false_L1_stress_to_L1_pressure_bound_used'])
    def test_kernel_near_far_coefficients_exact(self):self.assertEqual(kernel(),dict(near='-2/3',far='-2',integral='4',norm_decay='3/4'))
    def test_kernel_critical_endpoints_rejected(self):
        for p in (1,F(3,2),2):
            with self.assertRaises(ValueError):kernel(p)
    def test_cutoff_eight_retained_and_seven_fails(self):self.assertTrue(old.cutoff_weights(8)['valid']);self.assertFalse(old.cutoff_weights(7)['valid'])
    def test_quadratic_pressure_has_no_linear_B_half_term(self):self.assertEqual(len(flux_rows()['rows']),2);self.assertNotIn('unsplit_linear_pressure',[r['term'] for r in flux_rows()['rows']])
    def test_exact_Young_decays_four_and_fourteen_fifths(self):self.assertEqual(young_decay(F(3,2),1),4);self.assertEqual(young_decay(F(3,4),F(7,4)),F(14,5))
    def test_invalid_Young_exponents_rejected(self):
        for a,b in ((0,1),(2,1),(1,0)):
            with self.assertRaises(ValueError):young_decay(a,b)
    def test_new_full_radius_budget_two_and_old_four_thirds(self):self.assertEqual(flux_rows()['worst'],'2');self.assertEqual(flux_rows(False)['worst'],'4/3')
    def test_laplacian_scaling_limits_only_the_coarse_uniform_L2_estimate(self):self.assertTrue(exterior_review()['checks']['coarse_laplacian_flux_has_exact_radius_minus_two']);self.assertFalse(FLAGS['radius_minus_two_claimed_universal_optimal_rate_for_each_fixed_solution'])
    def test_fixed_T_precedes_radius_and_terminal_limits(self):self.assertIn('fix T<1, then R->infinity',protocol()['limits']);self.assertFalse(FLAGS['comparison_constants_uniform_through_t_one_assumed'])
    def test_zero_strain_needs_no_division(self):self.assertTrue(limit_review()['checks']['zero_strain_requires_no_division'])
    def test_protected_core_threshold_retained(self):self.assertIn('r0^2/(8*Xin)',limit_review()['core']);self.assertTrue(FLAGS['current_global_smooth_finite_energy_competitor_exclusion_relative_to_retained_lemmas'])
    def test_all_1015_inherited_and_18_direct_bindings(self):self.assertEqual(len(self.record['inherited_provenance']),1015);self.assertEqual(len(self.record['source_hashes']),18)
    def test_frozen_record_and_roundtrip_match(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r);validate(json.loads(json.dumps(self.record)))
    def test_unconditional_scope_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_cached_reviews_not_mutable_aliases(self):
        r=generate();r['reviews']['exact_complete_difference_energy']['checks']['complete_local_energy_flux_divergence_identity']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(local_review()['checks']['complete_local_energy_flux_divergence_identity'])
    def test_same_count_historical_key_substitution_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])
    def test_protocol_gates_sources_and_remaining_obligations_cannot_be_erased(self):
        for k in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[k],dict):r[k].pop(next(iter(r[k])))
            else:r[k].pop()
            with self.assertRaises(ValueError):validate(r)
    def test_cli_refuses_evidence_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')

if __name__=='__main__':unittest.main()
