#!/usr/bin/env python3
"""Current hierarchy, diagonal quantifiers, full cutoff controls and scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from realization import (B0,GAIN,KAPPA,RHO,stage_order,increment_orders,order_review,
    gain,tuple_loss,cutoff_loss,field_loss,comparison_budget,identity_review,scale_review)
from audit import HERE,FLAGS,generate,validate,compare_records,main,parent_review


class CurrentRealizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.record=generate()

    def test_current_seed_not_old_hypothetical_hierarchy(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'d065d09fa2f16710ea590454ebe5389e5edfb909')
        self.assertTrue(self.record['gates']['current_actual_finite_hierarchy_replayed'])

    def test_historical_realization_is_reproduced_not_assumed_current(self):
        r=self.record['reviews']['historical_realization_replayed_not_assumed_current']
        self.assertTrue(all(r['checks'].values()))
        self.assertTrue(r['old_conditional_hierarchy_not_substituted_for_current'])

    def test_previous_cycle_index_exact(self):
        self.assertEqual(B0-GAIN,RHO)
        self.assertEqual(increment_orders(1)['wave_velocity'],B0-KAPPA)
        self.assertNotEqual(increment_orders(1)['wave_velocity'],stage_order(1)-KAPPA)

    def test_affine_margin_holds_at_arbitrary_finite_stage(self):
        margins=order_review()['all_stage_margins_over_delta_j']
        for j in (1,2,1001,10**30):
            self.assertEqual({k:str(v-GAIN*j) for k,v in increment_orders(j).items()},margins)

    def test_complete_pressure_and_both_mean_representatives_present(self):
        rows=increment_orders(1)
        self.assertEqual(set(rows),{'wave_velocity','wave_potential','mean_azimuthal_potential',
            'direct_axisymmetric_swirl','complete_pressure_difference','mean_radial_velocity'})
        self.assertEqual(rows['complete_pressure_difference'],B0+Q(1,2)-2*KAPPA)

    def test_radial_mean_gain_and_smallest_margin(self):
        rows=increment_orders(17)
        self.assertEqual(rows['mean_radial_velocity']-rows['mean_azimuthal_potential'],1)
        self.assertEqual(min(v-GAIN*17 for v in rows.values()),1-4*KAPPA)

    def test_gain_is_current_cycle_gain_not_old_tenth(self):
        self.assertEqual(gain(7),Q(1,200)*7*(Q(1,2)-3*KAPPA))
        self.assertGreater(gain(7),Q(7,2000))

    def test_one_curl_and_fixed_log_absorption_reserved(self):
        for m in (0,2,30):
            self.assertEqual(cutoff_loss(m),tuple_loss(m+1))
            self.assertEqual(field_loss(m),cutoff_loss(m)+1)

    def test_invalid_indices_and_h_rejected(self):
        for j in (0,-1,True,1.5):
            with self.assertRaises(ValueError):increment_orders(j)
        for kwargs in ({'m':-1},{'m':True},{'target':-1},{'target':1.5},{'h':0},{'h':Q(1,100)}):
            with self.assertRaises(ValueError):comparison_budget(**kwargs)

    def test_fixed_stage_covers_all_residual_terms_and_jet_orders(self):
        r=self.record['reviews']['fixed_finite_stage_complete_residual_comparison']
        self.assertTrue(all(r['checks'].values()))
        self.assertEqual(len(r['rows']),12)

    def test_sample_comparison_is_only_a_rational_existence_budget(self):
        r=comparison_budget()
        self.assertEqual(r['comparison_stage'],12816)
        self.assertEqual(r['finite_stage_requirement'],1402)
        self.assertFalse(r['actual_fields_or_scales_evaluated'])

    def test_previous_stage_fails_sample_sufficient_tail_inequality(self):
        r=comparison_budget();J=r['comparison_stage']
        self.assertLess(gain(J)/2-cutoff_loss(2)-field_loss(2),1)
        self.assertGreaterEqual(gain(J+1)/2-cutoff_loss(2)-field_loss(2),1)

    def test_omitting_old_field_cross_is_not_a_valid_comparison(self):
        self.assertLess(Q(1,10)-field_loss(2),0)
        self.assertLess(gain(2)/2-cutoff_loss(2)-field_loss(2),0)

    def test_whole_interval_extreme_constant_recipe_uses_current_gain(self):
        r=scale_review();self.assertTrue(all(r['checks'].values()))
        self.assertEqual(r['rows'][0]['gamma'],str(gain(1)/2))
        self.assertGreater(r['rows'][-1]['log_scale'],10**14)
        self.assertTrue(r['synthetic_constants_only'])

    def controls(self,name):return identity_review()['controls'][name]['checks']

    def test_actual_implicit_q_in_all_third_order_jets(self):
        self.assertTrue(all(self.controls('actual_implicit_q').values()))

    def test_full_curl_product_rule_and_exact_divergence_controls(self):
        c=self.controls('complete_curl')
        self.assertTrue(c['full_curl_product_rule_through_two_velocity_derivatives'])
        self.assertTrue(c['complete_cut_velocity_is_divergence_free'])
        self.assertTrue(c['naive_velocity_cutoff_breaks_divergence'])

    def test_mean_radial_commutator_and_wrong_sign_controls(self):
        c=self.controls('complete_curl')
        self.assertTrue(c['mean_cutoff_adds_negative_axial_derivative_radial_velocity'])
        self.assertTrue(c['wrong_commutator_sign_is_detected'])

    def test_closed_core_axis_exterior_and_plateau_support_controls(self):
        c=self.controls('complete_curl')
        self.assertTrue(c['both_plateaus_are_exact'])
        self.assertTrue(c['core_axis_and_exterior_support_preserved'])

    def test_complete_residual_and_independent_cartesian_controls(self):
        c=self.controls('complete_residual')
        self.assertTrue(c['exact_complete_residual_difference'])
        self.assertTrue(c['independent_cartesian_residual_agrees'])
        self.assertTrue(c['every_linear_and_nonlinear_term_is_nonzero'])

    def test_cutoff_time_pressure_and_viscosity_cannot_be_dropped(self):
        c=self.controls('complete_residual')
        for n in ('cutoff_time_terms_cannot_be_dropped','cutoff_viscosity_commutators_cannot_be_dropped',
            'pressure_cutoff_gradient_cannot_be_dropped'):self.assertTrue(c[n])

    def test_raw_series_convergence_is_not_inferred_from_orders(self):
        self.assertTrue(all(self.controls('raw_series_control').values()))
        self.assertFalse(FLAGS['uniform_stage_constants_raw_uncut_series_or_fixed_q_contraction_verified'])

    def test_no_duplicate_roundoff_ledger_under_wrong_ancestor_comparator(self):
        r=identity_review()
        self.assertTrue(all(set(v)=={'checks'} for v in r['controls'].values()))

    def test_exact_Gaussian_tails_only_compared_at_fixed_finite_stage(self):
        self.assertTrue(FLAGS['designated_Gaussian_tails_never_inverted_or_independently_summed'])
        self.assertTrue(parent_review()[0]['flags']['actual_flat_velocities_and_exact_excluded_Gaussian_tail_ledger_retained'])

    def test_direct_swirl_not_promoted_to_compact_full_potential(self):
        self.assertFalse(FLAGS['direct_swirl_has_a_proved_compact_full_vector_potential'])
        self.assertTrue(order_review()['direct_swirl_is_not_assumed_to_have_a_compact_full_vector_potential'])

    def test_same_scale_sequence_not_j_of_q_or_uniform_constants(self):
        c=self.record['reviews']['exact_tail_pressure_and_quantifier_scope']['checks']
        for n in ('cutoff_scales_one_sequence_for_all_fixed_jet_orders_and_targets',
            'fixed_J_selected_before_q_limit_not_j_of_q',
            'no_uniform_stage_constant_or_raw_series_convergence_assumption'):self.assertTrue(c[n])

    def test_701_inherited_hashes_and_frozen_own_source_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),701)
        self.assertTrue(all(self.record['inherited_provenance'].values()))
        self.assertEqual(len(self.record['source_hashes']),15)

    def test_frozen_record_reproduces(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r)

    def test_record_roundtrips_without_numpy_boolean_serialization(self):
        r=json.loads(json.dumps(self.record));validate(r);compare_records(self.record,r)

    def test_unconditional_endpoint_global_force_and_blowup_promotions_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError):validate(r)

    def test_inherited_premises_cannot_be_erased(self):
        r=deepcopy(self.record)
        r['flags']['relies_on_inherited_leading_background_geometry_and_inverse_analytical_premises']=False
        with self.assertRaises(ValueError):validate(r)

    def test_changed_current_order_ledger_rejected_despite_true_gates(self):
        r=deepcopy(self.record)
        r['reviews']['current_all_stage_representatives_and_physical_losses']['delta']='1'
        with self.assertRaises(ValueError):validate(r)

    def test_changed_residual_budget_or_missing_curl_control_rejected(self):
        for name in ('budget','curl'):
            r=deepcopy(self.record)
            if name=='budget':r['reviews']['fixed_finite_stage_complete_residual_comparison']['rows'][0]['comparison_stage']=2
            else:r['reviews']['complete_implicit_q_curl_and_residual_controls']['controls']['complete_curl']['checks'].pop('naive_velocity_cutoff_breaks_divergence')
            with self.assertRaises(ValueError):validate(r)

    def test_missing_gate_obligation_or_provenance_rejected(self):
        for section in ('gates','remaining_obligations','source_hashes','inherited_provenance'):
            r=deepcopy(self.record)
            if section=='remaining_obligations':r[section].pop()
            else:r[section].pop(next(iter(r[section])))
            with self.assertRaises(ValueError):validate(r)

    def test_cli_never_overwrites_frozen_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':unittest.main()
