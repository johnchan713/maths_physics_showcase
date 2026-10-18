#!/usr/bin/env python3
"""Exact canonical force, all-support Borel and protected-plane controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
import audit
from force_bounds import universal_bounds,universal_checks,first_controlled,choose_dyadic_exponent,scale_review,protected_disk_review
from force_checks import protected_plane_review

class CanonicalForceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate()
    def test_current_canonical_endpoint_parent(self):self.assertEqual(audit.protocol()['parent_commit'],'99b5b2a72c16324a0c9389defa5e1d7a9cdfe76f')
    def test_canonical_input_all_gates(self):self.assertTrue(all(self.record['reviews']['canonical_current_field_input']['checks'].values()))
    def test_historical_force_not_identified_with_canonical(self):self.assertFalse(audit.FLAGS['historical_force_or_global_field_identified_with_canonical_field'])
    def test_twelve_full_operator_terms_retained(self):self.assertEqual(len(universal_bounds().TERMS),12)
    def test_pressure_cutoff_is_a_full_term(self):self.assertIn('pressure_cutoff',universal_bounds().TERMS)
    def test_force_Cauchy_orders_reserve_extra_time_derivative(self):
        for m in (0,1,5,100):
            r=universal_bounds().requirements(m);self.assertEqual(r['next_time_potential'],m+4);self.assertEqual(r['next_time_pressure'],m+2)
    def test_exact_localization_all_nine_vector_components(self):self.assertEqual(universal_checks().localization_review()['exact_localization_vector_components'],9)
    def test_all_twelve_omissions_detected(self):self.assertEqual(len(universal_checks().localization_review()['all_nonzero_omission_terms']),12)
    def test_naive_cutoff_breaks_divergence(self):self.assertTrue(universal_checks().localization_review()['checks']['naive_velocity_cutoff_breaks_divergence'])
    def test_nonaxisymmetric_direct_swirl_cutoff_breaks_divergence(self):self.assertTrue(universal_checks().localization_review()['checks']['theta_dependent_direct_swirl_cutoff_breaks_divergence'])
    def test_curl_commutator_orthogonality_is_exact(self):self.assertTrue(universal_checks().localization_review()['checks']['orthogonal_curl_commutator_has_zero_d_dot_grad_c'])
    def test_forced_pressure_Poisson_includes_nonzero_divergence(self):
        r=universal_checks().localization_review()['checks'];self.assertTrue(r['full_forced_pressure_Poisson_identity_exact']);self.assertTrue(r['omitted_force_divergence_is_nonzero'])
    def test_plane_fixture_all_288_mixed_component_controls(self):self.assertEqual(protected_plane_review()['exact_plane_component_comparisons'],288)
    def test_plane_fixture_not_promoted_to_infinite_flatness(self):self.assertTrue(protected_plane_review()['eighth_normal_jet_is_nonzero']);self.assertEqual(protected_plane_review()['normal_orders_through'],7)
    def test_full_analytical_disk_scope_not_only_origin(self):self.assertIn('equatorial disk',self.record['reviews']['protected_disk_and_future_force_flatness']['traces'])
    def test_future_uniform_normal_power_scope(self):self.assertIn('|z|^N',protected_disk_review()['future'])
    def test_flat_disk_not_an_open_force_free_neighborhood(self):self.assertFalse(audit.FLAGS['flat_disk_claimed_open_force_free_neighborhood'])
    def test_exact_integer_scale_acceptance_for_all_synthetic_norms(self):self.assertTrue(all(scale_review()['checks'].values()))
    def test_extreme_log2_norm_does_not_materialize_exponential(self):
        n=choose_dyadic_exponent(37,[(0,18,F(10**100,3))]);self.assertIsInstance(n,int);self.assertLessEqual(F(10**100,3)+(18-37)*n+37,-2)
    def test_scale_strictly_increases_from_previous_exponent(self):self.assertGreater(choose_dyadic_exponent(2,[(0,0,0)],1000),1001)
    def test_empty_constraints_rejected(self):
        with self.assertRaises(ValueError):choose_dyadic_exponent(1,[])
    def test_invalid_stage_and_previous_rejected(self):
        for j,previous in ((0,1),(-1,1),(True,1),(1,0),(1,True)):
            with self.assertRaises(ValueError):choose_dyadic_exponent(j,[(0,0,0)],previous)
    def test_negative_or_excess_derivative_constraints_rejected(self):
        for row in ((-1,0,0),(0,-1,0),(0,0,-1),(1,1,0),(0,2,0)):
            with self.assertRaises(ValueError):choose_dyadic_exponent(2,[row])
    def test_all_fixed_mixed_jets_eventually_controlled(self):
        for a,m in ((0,0),(1,3),(20,4),(6,70)):
            j=first_controlled(a,m);self.assertLessEqual(a+m,j//2);self.assertLess(m,j)
    def test_leibniz_all_98_polynomial_identities(self):self.assertEqual(universal_checks().Leibniz_review()['exact_Leibniz_identities'],98)
    def test_all_49_plateau_terminal_factorials(self):self.assertEqual(universal_checks().Leibniz_review()['exact_terminal_factorial_controls'],49)
    def test_no_full_force_Gevrey_or_Taylor_convergence_assumption(self):self.assertFalse(audit.FLAGS['complete_force_Gevrey_bound_or_analytic_Taylor_convergence_assumed'])
    def test_future_force_is_not_future_velocity(self):self.assertFalse(audit.FLAGS['future_velocity_or_pressure_constructed'])
    def test_pressure_and_transport_energy_work_zero(self):self.assertTrue(universal_checks().energy_review()['checks']['pressure_and_transport_energy_fluxes_exactly_zero'])
    def test_positive_viscous_dissipation_and_exact_energy_coefficients(self):self.assertTrue(universal_checks().energy_review()['checks']['viscous_work_equals_positive_dissipation']);self.assertTrue(universal_checks().energy_review()['checks']['integrated_energy_identity_has_factor_two_dissipation'])
    def test_initial_velocity_is_zero(self):self.assertTrue(universal_checks().energy_review()['checks']['manufactured_initial_velocity_exactly_zero'])
    def test_wrong_combined_extrema_budget_still_rejected(self):self.assertTrue(universal_bounds().energy_review()['checks']['separate_sup_and_total_bounds_not_added_with_same_budget'])
    def test_force_not_assumed_solenoidal(self):self.assertFalse(audit.FLAGS['force_assumed_divergence_free_or_compact_Leray_projection_used'])
    def test_current_official_source_date_and_force_requirements_retained(self):self.assertIn('checked 2026-10-18',audit.protocol()['official_source'])
    def test_frozen_record_exact_roundtrip(self):
        f=json.loads((audit.HERE/'evidence.json').read_text());audit.validate(f);audit.compare_records(self.record,f);audit.compare_records(self.record,json.loads(json.dumps(self.record)))
    def test_no_unconditional_or_future_velocity_scope_promotion(self):
        for flag in ('unconditional_admissible_force_residual_or_blowup_proved','future_velocity_or_pressure_constructed','canonical_global_competitor_pressure_recovery_and_comparison_joined','retained_core_background_and_finite_analytical_premises_discharged'):
            bad=deepcopy(self.record);bad['flags'][flag]=True
            with self.assertRaises(ValueError):audit.validate(bad)
    def test_source_hash_tamper_rejected(self):
        bad=deepcopy(self.record);key=next(iter(bad['source_hashes']));bad['source_hashes'][key]='0'*64
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_same_count_provenance_substitution_rejected(self):
        bad=deepcopy(self.record);v=bad['inherited_provenance']['bindings'];v['forged']=v.pop(next(iter(v)))
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_exact_scale_margin_tampering_rejected(self):
        bad=deepcopy(self.record);bad['reviews']['exact_rational_Borel_scale']['rows'][0]['log2_margin_upper']='1'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_cached_plane_reviews_not_mutated(self):
        bad=deepcopy(self.record);bad['reviews']['finite_order_protected_plane_controls']['checks'].clear()
        self.assertTrue(audit.reviews()['finite_order_protected_plane_controls']['checks'])
    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'kept.json';p.write_text('keep')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
if __name__=='__main__':unittest.main()
