#!/usr/bin/env python3
"""Exact complete localization, force jets and compact energy controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import tempfile
import unittest
from localized_bounds import TERMS,order,requirements,first_controlled,derivative_factor,domain_review,terminal_review,extension_review,energy_review,scope_review
from localized_checks import localization_review,Leibniz_review,energy_review as exact_energy
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main

class ForceFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_endpoint_parent_is_exact(self):self.assertEqual(protocol()['parent_commit'],'4899fa23a49f654c53bd3332c4c269bf9072e744')
    def test_current_enlarged_margin_not_historical_narrow_margin(self):self.assertIn('16q<=6q_e',protocol()['enlarged_observation']);self.assertLess(6*F(1,16),1)
    def test_localization_keeps_curl_commutator_and_full_pressure(self):self.assertIn('grad(c) cross A',protocol()['localization']);self.assertIn('full fixed gauge',protocol()['localization'])
    def test_twelve_named_localization_terms_are_distinct(self):self.assertEqual(len(TERMS),12);self.assertEqual(len(set(TERMS)),12)
    def test_force_jet_orders_reserve_three_potential_and_one_pressure(self):self.assertEqual(requirements(5),dict(force=5,potential=8,direct_swirl=7,pressure=6,next_time_potential=9,next_time_direct_swirl=8,next_time_pressure=7))
    def test_invalid_finite_orders_rejected(self):
        for n in (-1,True,F(1,2)):
            with self.assertRaises(ValueError):order(n)
            with self.assertRaises(ValueError):requirements(n)
    def test_exact_localization_at_three_positive_viscosities(self):self.assertEqual(localization_review()['exact_localization_vector_components'],9);self.assertTrue(localization_review()['checks']['all_twelve_localization_terms_exact'])
    def test_every_one_of_twelve_omissions_is_nonzero(self):self.assertEqual(len(localization_review()['all_nonzero_omission_terms']),12)
    def test_complete_cut_potential_is_divergence_free(self):self.assertTrue(localization_review()['checks']['complete_cut_potential_velocity_is_divergence_free'])
    def test_naive_velocity_cutoff_fails(self):self.assertTrue(localization_review()['checks']['naive_velocity_cutoff_breaks_divergence'])
    def test_nonaxisymmetric_direct_swirl_cutoff_fails(self):self.assertTrue(localization_review()['checks']['theta_dependent_direct_swirl_cutoff_breaks_divergence'])
    def test_curl_commutator_orthogonal_to_cutoff_gradient(self):self.assertTrue(localization_review()['checks']['orthogonal_curl_commutator_has_zero_d_dot_grad_c'])
    def test_forced_pressure_Poisson_keeps_force_divergence(self):self.assertTrue(localization_review()['checks']['full_forced_pressure_Poisson_identity_exact']);self.assertTrue(localization_review()['checks']['omitted_force_divergence_is_nonzero'])
    def test_all_terminal_force_jets_flat_at_origin_not_velocity(self):self.assertTrue(terminal_review()['checks']['no_regular_velocity_or_pressure_at_origin_inferred']);self.assertIn('every spatial jet at origin zero',terminal_review()['trace'])
    def test_all_fixed_jet_pairs_enter_one_sequence(self):
        for a,m in ((0,0),(1,3),(20,4),(6,70)):
            j=first_controlled(a,m);self.assertLessEqual(a+m,j//2);self.assertLess(m,j)
    def test_Leibniz_factor_keeps_factorial(self):self.assertEqual(derivative_factor(5,2,1),F(1,12));self.assertEqual(derivative_factor(3,2,1),1)
    def test_invalid_Leibniz_terms_are_zero(self):self.assertEqual(derivative_factor(1,3,1),0);self.assertEqual(derivative_factor(3,1,2),0)
    def test_98_exact_whole_transition_Leibniz_identities(self):self.assertEqual(Leibniz_review()['exact_Leibniz_identities'],98)
    def test_49_exact_plateau_terminal_controls(self):self.assertEqual(Leibniz_review()['exact_terminal_factorial_controls'],49);self.assertTrue(Leibniz_review()['checks']['missing_factorial_changes_terminal_jet'])
    def test_whole_support_Borel_exponent_not_just_transition(self):self.assertIn('m-j',extension_review()['derivative_bound']);self.assertTrue(Leibniz_review()['checks']['uniform_whole_support_exponent_is_m_minus_j'])
    def test_extreme_norm_constraints_and_increasing_scales(self):
        r=extension_review();self.assertTrue(all(x['log_margin']<-1 and x['increased'] for x in r['rows']));self.assertTrue(r['actual_terminal_norms_and_extension_scales_not_evaluated'])
    def test_geometric_tail_constant_has_factor_two(self):self.assertIn('2^(1-J)',extension_review()['tail'])
    def test_no_future_velocity_from_smooth_force(self):self.assertFalse(FLAGS['future_global_velocity_or_pressure_constructed']);self.assertFalse(FLAGS['analytic_or_zero_heat_extension_used'])
    def test_compact_force_support_and_zero_initial_datum(self):self.assertEqual(scope_review()['data'],'u(0)=0');self.assertIn('[1-tau0,3/2]',scope_review()['force'])
    def test_no_unforced_pressure_projection(self):self.assertFalse(FLAGS['force_assumed_divergence_free']);self.assertIn('-div f',protocol()['pressure'])
    def test_exact_compact_pressure_and_transport_energy_fluxes(self):self.assertTrue(exact_energy()['checks']['pressure_and_transport_energy_fluxes_exactly_zero'])
    def test_exact_differential_and_integrated_energy_factors(self):
        r=exact_energy()['checks'];self.assertTrue(r['differential_energy_identity_with_half_coefficient_exact']);self.assertTrue(r['integrated_energy_identity_has_factor_two_dissipation'])
    def test_exact_viscous_work_positive_and_initial_zero(self):
        r=exact_energy();self.assertGreater(F(r['exact_dissipation_coefficient']),0);self.assertTrue(r['checks']['manufactured_initial_velocity_exactly_zero'])
    def test_zero_energy_start_uses_regularized_square_root(self):self.assertIn('delta',energy_review()['regularization'])
    def test_separate_extrema_cannot_share_one_energy_budget(self):
        c=energy_review()['wrong_combined_budget_counterexample'];self.assertGreater(F(c['lower_bound_for_sup_energy_plus_twice_total_dissipation']),F(c['force_budget_squared']))
    def test_finite_regularity_polynomials_not_actual_C_infinity_fields(self):self.assertTrue(exact_energy()['finite_regularity_fixture_not_the_actual_C_infinity_field']);self.assertTrue(localization_review()['manufactured_local_polynomial_fields_not_actual_compact_C_infinity_flow'])
    def test_all_997_historical_and_18_direct_bindings(self):self.assertEqual(len(self.record['inherited_provenance']),997);self.assertEqual(len(self.record['source_hashes']),18)
    def test_frozen_record_and_roundtrip_match(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r);validate(json.loads(json.dumps(self.record)))
    def test_scope_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)
    def test_cached_reviews_not_mutable_aliases(self):
        r=generate();r['reviews']['exact_complete_localization_and_pressure']['checks']['all_twelve_localization_terms_exact']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(localization_review()['checks']['all_twelve_localization_terms_exact'])
    def test_same_count_historical_key_substitution_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])
    def test_protocol_gates_sources_and_obligations_cannot_be_erased(self):
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
