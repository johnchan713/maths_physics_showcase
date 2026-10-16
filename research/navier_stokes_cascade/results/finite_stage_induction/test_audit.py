#!/usr/bin/env python3
"""All-half-line closure, sharp gain controls, current seed and quantifiers."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,generate,validate,compare_records,main,parent_review
from induction import (Affine,line,certificate,seed_review,stage_order,physical_loss,
    physical_power,first_stage_for_power,harmonic_bound,B0,RHO,MU,GAIN,KAPPA)


class FiniteInductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.record=generate()

    def review(self,name): return self.record['reviews'][name]

    def test_current_seed_is_latest_actual_state(self):
        r=self.review('current_actual_seed')
        self.assertEqual(Q(r['seed_orders']['B']),B0)
        self.assertEqual(Q(r['seed_orders']['C']),B0+Q(1,2))
        self.assertTrue(r['checks']['actual_latest_seed_not_old_reference_initialization'])

    def test_current_defects_stronger_than_induction_input(self):
        self.assertGreater(Q(self.review('current_actual_seed')['seed_margins']['defects']),0)

    def test_old_reference_seed_not_silently_substituted(self):
        r=deepcopy(parent_review()[0]);r['reviews']['current_mean_gain_physical_and_uncancelled_wave_ledger']['rows']['final_nonzero_regular']='7/10'
        self.assertFalse(all(seed_review(r)['checks'].values()))

    def test_current_seed_requires_actual_tail_exclusion(self):
        r=deepcopy(parent_review()[0]);r['flags']['actual_flat_velocities_in_every_product_and_Gaussian_tail_exclusion_preserved']=False
        self.assertFalse(seed_review(r)['checks']['exact_Gaussian_tail_exclusion_is_a_seed_input'])

    def test_affine_boundary_margin_is_exact(self):
        self.assertEqual(line(2,-KAPPA).margin(line(1,GAIN)),B0-KAPPA-GAIN)

    def test_negative_slope_cannot_be_certified_by_positive_boundary_value(self):
        with self.assertRaises(ValueError): line(0,1000).margin(line(1,0))

    def test_all_half_line_margins_nonnegative(self):
        c=certificate();self.assertTrue(all(v>=0 for v in c['half_line_margins'].values()))

    def test_three_saturated_rows_identify_sharp_available_gain(self):
        self.assertEqual(set(certificate()['limiting_rows']),{
            'nonzero_regular:sourced_full_linear','nonzero_regular:signed_full_linear','averaged_signed:old_complex_error_cross'})

    def test_arbitrarily_small_gain_increase_is_rejected_exactly(self):
        with self.assertRaisesRegex(ValueError,'sourced_full_linear.*old_complex_error_cross'):
            certificate(GAIN+Q(1,10**80))

    def test_nonpositive_gain_rejected_and_smaller_positive_gain_valid(self):
        for g in (0,-1):
            with self.assertRaises(ValueError): certificate(g)
        self.assertGreater(certificate(GAIN/2)['minimum_margin'],0)

    def test_wave_linear_viscous_and_pressure_remainder_row(self):
        c=certificate();self.assertEqual(c['groups']['nonzero_regular']['sourced_full_linear'],line(1,Q(1,2)-3*KAPPA))

    def test_old_new_and_new_square_wave_terms_all_present(self):
        rows=certificate()['groups']['nonzero_regular']
        for n in ('sourced_old_wave_cross','sourced_new_square','sourced_actual_mean_cross',
            'signed_old_wave_cross','signed_new_square','signed_actual_mean_cross','fast_mean_wave_change','slow_mean_wave_change'):
            self.assertIn(n,rows)

    def test_generic_old_complex_cross_limits_mean_target(self):
        c=certificate();self.assertEqual(c['groups']['averaged_signed']['old_complex_error_cross'],line(1,RHO))
        self.assertEqual(RHO-Q(1,2),GAIN)

    def test_full_new_curl_cross_signed_square_and_axial_pressure_retained(self):
        rows=certificate()['groups']['averaged_signed']
        self.assertEqual(set(rows),{'old_complex_error_cross','complete_old_new_curl','signed_square','retained_current_moment_bump','full_axial_flux_and_pressure'})

    def test_fast_viscosity_exceeds_final_mean_target_by_half(self):
        c=certificate();self.assertEqual(c['half_line_margins']['fast_full_mean:radial_viscosity'],Q(1,2))

    def test_accumulated_old_Haar_mean_limits_slow_scalar_gain(self):
        c=certificate();rows=c['groups']['postslow_defects']
        self.assertEqual(min(v.at(B0) for v in rows.values()),B0+Q(3,2))
        self.assertLess(B0+Q(3,2),2*(B0+Q(1,2)))

    def test_fixed_wave_and_total_mean_orders_preserved(self):
        c=certificate();self.assertTrue(all(v>0 for v in c['invariant_margins'].values()))
        self.assertGreater(B0,RHO);self.assertGreater(B0+Q(1,2)-KAPPA,MU)

    def test_current_postfast_defects_remain_valid_slow_inputs(self):
        self.assertEqual(certificate()['invariant_margins']['postfast_defect_change_over_C'],MU-KAPPA)

    def test_exact_stage_formula_for_arbitrary_finite_integer(self):
        self.assertEqual(stage_order(10**30),B0+GAIN*10**30)
        self.assertEqual(stage_order(1),Q(2)-9*KAPPA)

    def test_stage_validation_rejects_nonintegers_and_boolean(self):
        for j in (-1,True,1.0,Q(1,2)):
            with self.assertRaises(ValueError): stage_order(j)

    def test_finite_harmonic_recurrence_without_enumerating_modes(self):
        self.assertEqual(harmonic_bound(20),4*2**20)
        self.assertEqual(harmonic_bound(21),2*harmonic_bound(20))

    def test_physical_loss_is_independent_of_stage(self):
        self.assertEqual(physical_loss(Q(1,200),3),Q(851,100))
        self.assertEqual(physical_power(100,Q(1,200),3)-physical_power(99,Q(1,200),3),Q(1,200)*GAIN)

    def test_physical_parameter_and_jet_validation(self):
        for h,m in ((0,0),(Q(1,100),0),(Q(1,200),-1),(Q(1,200),True)):
            with self.assertRaises(ValueError): physical_loss(h,m)

    def test_low_stage_coarse_physical_power_remains_negative(self):
        self.assertEqual(physical_power(0,Q(1,200),0),Q(-25025003,10000000))

    def test_finite_target_stage_threshold_is_minimal(self):
        j=first_stage_for_power(0)
        self.assertEqual(j,1002)
        self.assertGreaterEqual(physical_power(j,Q(1,200),0),0)
        self.assertLess(physical_power(j-1,Q(1,200),0),0)
        self.assertEqual(first_stage_for_power(-100),0)

    def test_each_fixed_derivative_target_has_finite_order_stage(self):
        for target,m in ((1,2),(10,5),(100,20)):
            j=first_stage_for_power(target,m=m)
            self.assertGreaterEqual(physical_power(j,Q(1,200),m),target)
            self.assertLess(physical_power(j-1,Q(1,200),m),target)

    def test_arbitrary_stage_constants_defeat_a_fixed_q_contraction_claim(self):
        r=self.review('stage_independent_physical_loss_and_quantifier_boundary')
        self.assertTrue(r['checks']['permitted_constant_growth_can_defeat_fixed_q_bound_contraction'])
        self.assertIn('not actual residual data',r['admissible_constant_growth_example'])

    def test_quantifiers_do_not_promote_to_single_infinite_field(self):
        r=self.review('stage_independent_physical_loss_and_quantifier_boundary')
        self.assertTrue(r['no_fixed_q_contraction_uniform_stage_constants_or_j_of_q_substitution'])
        self.assertFalse(FLAGS['single_infinite_realized_field_with_flat_complete_residual_verified'])

    def test_complete_interaction_and_missing_hypothesis_controls_replayed(self):
        r=self.review('replayed_complete_interaction_and_missing_hypothesis_controls')
        self.assertTrue(all(r['checks'].values()));self.assertEqual(len(r['checks']),5)

    def test_parent_and_686_inherited_hashes(self):
        self.assertTrue(self.record['gates']['current_complex_mean_parent_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),686)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_evidence_reproduces(self):
        frozen=json.loads((HERE/'evidence.json').read_text());validate(frozen);compare_records(self.record,frozen)

    def test_unconditional_single_field_force_and_blowup_promotions_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError): validate(r)

    def test_conditional_premises_and_cycle_order_cannot_change(self):
        for section,n,v in (('protocol','cycle_order','forward invert Gaussian tails'),
            ('flags','relies_on_replayed_leading_background_geometry_and_inverse_analytical_premises',False)):
            r=deepcopy(self.record);r[section][n]=v
            with self.assertRaises(ValueError): validate(r)

    def test_missing_gate_obligation_or_provenance_rejected(self):
        for section in ('gates','remaining_obligations','source_hashes','inherited_provenance'):
            r=deepcopy(self.record)
            if section=='remaining_obligations':r[section].pop()
            else:r[section].pop(next(iter(r[section])))
            with self.assertRaises(ValueError):validate(r)

    def test_changed_affine_ledger_cannot_hide_behind_passed_gates(self):
        r=deepcopy(self.record)
        r['reviews']['complete_four_step_half_line_induction']['certificate']['gain']='1'
        with self.assertRaises(ValueError):validate(r)

    def test_cli_does_not_overwrite_frozen_record(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__': unittest.main()
