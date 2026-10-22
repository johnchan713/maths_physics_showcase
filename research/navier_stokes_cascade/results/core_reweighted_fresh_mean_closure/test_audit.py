#!/usr/bin/env python3
"""Same-gap higher-pole and post-sourced full-pressure mean regression checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review
from budgets import P,cut_recipe,shift_recipe,materialize_shift,translate


class FreshMeanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        rows=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(rows);self.assertTrue(all(rows))

    def test_fixed_shift_covers_every_new_operator_domain(self):self.control('one_fixed_shift_covers_all_four_new_domains')
    def test_same_gap_recipe_matches_existing_independent_materializer(self):self.control('same_actual_gap_kernel_recipe')
    def test_original_width_and_gap_are_preserved(self):self.control('same_original_width_and_gap')
    def test_exact_translation_works_below_shift_without_new_band_gate(self):self.control('exact_binomial_translation_no_new_S_ge_K_gate')
    def test_positive_polynomials_majorize_original_scale_inputs(self):self.control('positive_polynomial_inputs_monotone_under_estimate_shift')
    def test_giant_constant_recipe_does_not_need_materialization(self):self.control('giant_constant_retained_exactly_without_materialization')
    def test_diagnostic_exponent_guard(self):self.control('materializer_guard_is_explicit_diagnostic_limit')
    def test_old_slow_and_fast_means_are_accumulated(self):self.control('old_accumulated_means_survive_fresh_fast_and_slow')
    def test_changed_complete_wave_keeps_old_crosses_and_sourced_square(self):self.control('complete_changed_wave_keeps_old_cross_and_sourced_square')
    def test_fresh_postfast_full_pressure_and_scalars(self):self.control('all_three_scalars_and_pressure_recomputed_after_fast')
    def test_slow_uses_fresh_pressure_and_all_three_scalar_targets(self):self.control('slow_map_uses_fresh_postfast_pressure_and_all_three_scalars')
    def test_postslow_pressure_is_reconstructed_again(self):self.control('final_pressure_and_scalars_recomputed_after_slow')
    def test_postslow_pressure_includes_slow_products(self):self.control('postslow_pressure_contains_slow_mean_products')
    def test_full_pressure_and_moment_envelopes_are_reweighted(self):self.control('full_pressure_Kp_and_global_moment_Kmom_reweighted')
    def test_axial_A1_cost_is_kept(self):self.control('axial_A1_operation_retained_in_fast_tangential_bound')
    def test_finite_composed_jet_loss(self):self.control('postslow_order_zero_uses_finite_original_mean_jets_through_31')
    def test_generated_full_cut_pressure_is_present(self):self.control('full_cut_sourced_pressure_retained')
    def test_original_normal_and_tangent_budget_is_used(self):self.control('normal_and_tangent_from_original_frame_inputs')
    def test_wrong_next_source_harmonic_grammar_is_blocked(self):
        for name in ('source','raw_source','inverse_source'):self.control('next_harmonic_grammar_blocks_'+name)

    def test_invalid_poles_geometry_and_approximate_inputs_rejected(self):
        for values in ((-1,4,Q(1,2),1),(1,0,Q(1,2),1),(1,4,0,1),(1,4,2,1),(1,4,Q(1,2),0),(1,4,Q(1,2),3),(0.5,4,Q(1,2),1)):
            with self.assertRaises((ValueError,TypeError)):cut_recipe(*values)
        with self.assertRaises(ValueError):translate(P(1),-1)
        with self.assertRaises(TypeError):translate(1,1)
        for x in (-1,0.5):
            with self.assertRaises((ValueError,TypeError)):materialize_shift(shift_recipe(1,4,Q(1,2),1),x)

    def test_exact_shift_translation_keeps_each_mixed_power(self):
        f=translate(P({2:1}),Q(3,2))
        self.assertEqual(f.c,{0:Q(9,4),1:3,2:1})

    def test_four_records_no_new_assumptions_and_no_actual_claim(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),4)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertTrue(all(self.record['gates'].values()));self.control('no_actual_input_enclosure_or_residual_promotion')

    def test_frozen_record_reproduces(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_every_scope_flag_and_obligation_is_bound(self):
        for k,v in self.record['flags'].items():
            r=deepcopy(self.record);r['flags'][k]=not v
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)

    def test_frozen_output_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
