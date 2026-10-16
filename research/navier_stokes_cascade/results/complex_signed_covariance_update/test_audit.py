#!/usr/bin/env python3
"""Generic covariance, full pressure, restrictive replay repair and scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,generate,validate,compare_records,main
from covariance import KAPPA,orders,physical_powers
from replay import compatible_compare,ROUND_FIELDS,ABS_ROUNDOFF


class ComplexSignedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.record=generate()

    def review(self,name): return self.record['reviews'][name]
    def check(self,review,name): self.assertTrue(self.review(review)['checks'][name])

    def test_generic_full_tensor_exact_identity_and_independent_average(self):
        self.check('generic_complete_tensor_and_complex_old_error','exact_generic_full_tensor_increment')
        self.assertEqual(self.review('generic_complete_tensor_and_complex_old_error')['exact_tensor_checks'],36)

    def test_complete_old_and_new_curls_divergence(self):
        self.check('generic_complete_tensor_and_complex_old_error','old_and_new_complete_curls_are_divergence_free')
        self.assertEqual(self.review('generic_complete_tensor_and_complex_old_error')['exact_divergence_checks'],8)

    def test_each_generic_error_group_is_needed(self):
        self.check('generic_complete_tensor_and_complex_old_error','all_three_generic_remainder_groups_are_needed')

    def test_old_complex_error_prevents_initial_real_shortcut(self):
        self.check('generic_complete_tensor_and_complex_old_error','generic_remainder_not_the_initial_real_formula')

    def test_full_six_entry_tensor_not_only_cross_flux(self):
        for n in ('all_six_covariance_changes_retained','whole_tensor_identity'):
            self.check('generic_complete_tensor_and_complex_old_error',n)

    def test_compact_recomputed_pressure(self):
        self.check('recomputed_full_axial_pressure_and_defects','full_compact_pressure')

    def test_axial_pressure_moment_retains_whole_bump(self):
        self.check('recomputed_full_axial_pressure_and_defects','complete_axial_pressure_moment')

    def test_both_integrated_fluxes(self):
        for n in ('exact_angular_integrated_flux','exact_axial_integrated_flux'):
            self.check('recomputed_full_axial_pressure_and_defects',n)

    def test_target_cancellation_keeps_all_nonradial_pressure_and_flux_terms(self):
        for n in ('angular_target_cancellation','axial_target_cancellation'):
            self.check('recomputed_full_axial_pressure_and_defects',n)

    def test_axial_pressure_and_moving_bump_cannot_be_dropped(self):
        self.check('recomputed_full_axial_pressure_and_defects','axial_pressure_covariance_and_moving_bump_terms_nonzero')

    def test_all_three_defects_recomputed_and_nonzero(self):
        self.check('recomputed_full_axial_pressure_and_defects','all_three_updated_defects_are_nonzero')

    def test_actual_fixed_inverse_factor_two_and_signed_target(self):
        for n in ('fixed_signed_linear_inverse_recovers_target','one_inverse_coefficient_negative','factor_two_is_necessary'):
            self.check('current_fixed_signed_inverse',n)

    def test_squared_partition_normalization_without_dividing_by_partition(self):
        for n in ('squared_partition_physical_normalization','no_division_by_slow_partition_needed'):
            self.check('current_fixed_signed_inverse',n)

    def test_half_edge_weight_and_no_updated_positivity_requirement(self):
        for n in ('every_checked_reciprocal_root_term_retains_half_edge_weight','weighted_class_does_not_ensure_positive_updated_amplitudes'):
            self.check('current_fixed_signed_inverse',n)

    def test_current_Haar_target_not_unaveraged_graph_order(self):
        r=self.review('current_generic_orders_and_physical_sources')['rows']
        self.assertEqual(Q(r['target']),Q(3,2)-3*KAPPA)
        self.assertEqual(Q(r['old_remainder']),1-3*KAPPA)

    def test_old_error_cross_and_square_both_limit_generic_gain(self):
        r=self.review('current_generic_orders_and_physical_sources')['rows']
        self.assertEqual(Q(r['old_error_new_principal']),2-6*KAPPA)
        self.assertEqual(Q(r['signed_square']),2-6*KAPPA)
        self.assertGreater(Q(r['complete_old_new_curl']),Q(r['generic_error']))

    def test_generic_gain_condition_at_and_above_boundaries(self):
        alpha=Q(3,2)-3*KAPPA;rho=1-3*KAPPA
        self.assertEqual(orders(alpha,Q(1,2))['gain'],0)
        self.assertEqual(orders(1,rho)['gain'],0)
        for a,r in ((Q(11,10),Q(3,5)),(Q(3,2),1),(2,2)):
            self.assertEqual(orders(a,r)['gain'],min(r-Q(1,2),Q(1,2)-KAPPA,a-1))
            self.assertGreater(orders(a,r)['gain'],0)

    def test_invalid_generic_orders_rejected(self):
        for a,r,k in ((Q(9,10),1,KAPPA),(2,Q(2,5),KAPPA),(2,1,0),(2,1,Q(1,2))):
            with self.assertRaises(ValueError): orders(a,r,k)

    def test_full_tensor_change_does_not_gain_to_generic_pair_error(self):
        r=orders(Q(3,2)-3*KAPPA,1-3*KAPPA)
        self.assertEqual(r['full_covariance_change'],r['target'])
        self.assertLess(r['full_covariance_change'],r['generic_error'])

    def test_full_mean_and_defects_wait_for_relaxation(self):
        r=self.review('current_generic_orders_and_physical_sources')['rows']
        self.assertEqual(Q(r['full_mean_and_pressure']),Q(3,2)-4*KAPPA)
        self.assertEqual(Q(r['defects']),Q(3,2)-3*KAPPA)
        self.assertEqual(Q(r['new_averaged_tangential']),2-6*KAPPA)

    def test_nonzero_regular_residual_preserved_not_improved(self):
        r=self.review('current_generic_orders_and_physical_sources')
        self.assertEqual(min(map(Q,r['generated_nonzero_orders'].values())),Q(3,2)-6*KAPPA)
        self.assertEqual(r['rows']['old_and_new_nonzero_residual'],r['rows']['previous_full_wave_residual'])

    def test_current_primary_remainder_and_finite_harmonics(self):
        r=self.review('current_generic_orders_and_physical_sources')
        self.assertEqual(Q(r['rows']['current_fixed_primary_remainder']),1-3*KAPPA)
        self.assertEqual(r['maximum_updated_residual_multiplier'],4)

    def test_all_fixed_physical_derivative_powers(self):
        a=Q(3,2)-3*KAPPA;h=Q(1,200);p=physical_powers(a,h,0,0)
        self.assertEqual(p,dict(velocity=Q(-10000003,20000000),pressure=Q(-20050003,20000000),potential=Q(49997,20000000)))
        for x,t in ((1,0),(0,1),(3,2),(6,4)):
            shift=-x*Q(201,400)-t*Q(403,400)
            self.assertEqual(physical_powers(a,h,x,t),{n:v+shift for n,v in p.items()})

    def test_physical_order_and_derivative_validation(self):
        for a,h,x,t in ((0,Q(1,200),0,0),(2,0,0,0),(2,Q(1,100),0,0),(2,Q(1,200),-1,0),(2,Q(1,200),True,0)):
            with self.assertRaises(ValueError): physical_powers(a,h,x,t)

    def test_exactly_three_roundoff_fields_with_stricter_ceiling(self):
        self.assertEqual(len(ROUND_FIELDS),3);self.assertLess(ABS_ROUNDOFF,2e-10)
        self.check('parent_CI_roundoff_compatibility_not_scope_relaxation','adapter_stricter_than_unchanged_fixture_gates')

    def test_roundoff_change_accepted_only_at_listed_paths(self):
        strict=lambda a,b,path: self.assertEqual(a,b)
        for field in ROUND_FIELDS:
            compatible_compare(8e-15,5e-15,field,strict)
            compatible_compare(0.,5e-15,'root.reviews.'+field,strict)
        with self.assertRaises(AssertionError): compatible_compare(8e-15,5e-15,'unrelated.error',strict)

    def test_adapter_rejects_nonfinite_negative_and_excessive_error(self):
        field=next(iter(ROUND_FIELDS))
        for bad in (-1e-15,3e-12,float('nan'),float('inf')):
            with self.assertRaises(ValueError): compatible_compare(bad,5e-15,field,lambda *_:None)

    def test_adapter_does_not_relax_types_keys_or_list_lengths(self):
        for a,b in ((True,1),({'x':1},{}),([1],[])):
            with self.assertRaises(ValueError): compatible_compare(a,b,'root',lambda *_:None)

    def test_adapter_cannot_hide_false_scope_or_omitted_term(self):
        self.check('parent_CI_roundoff_compatibility_not_scope_relaxation','all_other_fields_keep_inherited_comparison')
        self.assertFalse(FLAGS['parent_CI_claimed_successful'])

    def test_diagnostics_are_not_actual_global_fields(self):
        r=self.review('generic_complete_tensor_and_complex_old_error')
        self.assertTrue(r['manufactured_rotated_old_curl_not_actual_sourced_pulse'])
        self.assertTrue(r['actual_closed_label_separation_not_replaced_by_fixture_frequencies'])

    def test_parent_and_655_inherited_hashes(self):
        self.assertTrue(self.record['gates']['actual_first_sourced_parent_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),655)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_evidence_reproduces(self):
        frozen=json.loads((HERE/'evidence.json').read_text());validate(frozen);compare_records(self.record,frozen)

    def test_false_scope_promotions_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError): validate(r)

    def test_changed_CI_history_or_analytical_dependencies_rejected(self):
        for section,n,v in (('protocol','parent_CI_conclusion','success'),
            ('flags','relies_on_replayed_background_geometry_pulse_and_mean_analytical_lemmas',False)):
            r=deepcopy(self.record);r[section][n]=v
            with self.assertRaises(ValueError): validate(r)

    def test_changed_ledger_missing_gate_obligation_or_provenance_rejected(self):
        r=deepcopy(self.record);r['reviews']['current_generic_orders_and_physical_sources']['rows']['generic_error']='3'
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
