#!/usr/bin/env python3
"""All-finite current-source half-line, interaction and quantifier controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from finite_bounds import KAPPA,B0,RHO,MU,GAIN,line,certificate,stage_order,harmonic_bound,next_harmonics,physical_loss,physical_power,first_stage
from finite_checks import harmonic_review,covariance_review
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main


class FiniteInductionFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_parent_and_unchanged_operations(self):
        self.assertEqual(protocol()['parent_commit'],'89f9c856d8800a71f16c36e23819babb61aa0c63')
        self.assertIn('operation order',protocol()['unchanged_parameters'])

    def test_current_seed_B0_full_mean_and_stronger_defects(self):
        r=self.record['reviews']['current_actual_seed'];self.assertTrue(all(r['checks'].values()))
        self.assertEqual(F(r['orders']['B']),B0);self.assertEqual(F(r['orders']['C']),B0+F(1,2))

    def test_exact_sharper_gain_rho_and_mu(self):
        self.assertEqual(GAIN,F(1,2)-KAPPA);self.assertEqual(RHO,MU);self.assertEqual(RHO-F(1,2),GAIN)

    def test_all_thirty_one_complete_rows_on_entire_half_line(self):
        c=certificate();self.assertEqual(len(c['margins']),31);self.assertGreaterEqual(min(c['margins'].values()),0)
        self.assertIn('every B>=',c['proof_domain'])

    def test_negative_affine_slope_rejected_even_if_boundary_positive(self):
        with self.assertRaises(ValueError):line(0,1000).margin(line(1,0))

    def test_boundary_margin_of_nondecreasing_difference_is_exact(self):
        self.assertEqual(line(2,1).margin(line(1,0)),B0+1)

    def test_seven_regular_rows_limit_gain_including_fast_mean_wave(self):
        rows=certificate()['limiting_rows'];self.assertEqual(len([n for n in rows if n.startswith('nonzero:')]),7)
        self.assertIn('nonzero:fast_mean_wave',rows)

    def test_both_generic_covariance_crosses_limit_mean_gain(self):
        rows=certificate()['limiting_rows'];self.assertEqual({n for n in rows if n.startswith('averaged_signed:')},{'averaged_signed:old_complex_cross','averaged_signed:complete_old_new_curl'})

    def test_larger_gain_rejected_by_available_estimates(self):
        with self.assertRaisesRegex(ValueError,'sourced_full_linear'):certificate(GAIN+KAPPA)

    def test_zero_or_negative_gain_rejected(self):
        for g in (0,-1):
            with self.assertRaises(ValueError):certificate(g)

    def test_conservative_old_gain_still_accepted(self):
        c=certificate(F(1,2)-3*KAPPA);self.assertGreater(min(c['margins'].values()),0)

    def test_fast_viscosity_and_postslow_defect_margins(self):
        c=certificate();self.assertEqual(c['margins']['fast_full_mean:radial_viscosity'],F(1,2)-2*KAPPA)
        self.assertEqual(c['margins']['postslow_defects:accumulated_Haar_cross'],F(1,2)+KAPPA)

    def test_all_preserved_state_invariants_strict(self):self.assertGreater(min(certificate()['invariants'].values()),0)

    def test_stage_formula_and_exact_positive_increment(self):
        self.assertEqual(stage_order(0),B0);self.assertEqual(stage_order(23)-stage_order(22),GAIN)

    def test_invalid_stage_indices_rejected(self):
        for j in (-1,True,F(1,2)):
            with self.assertRaises(ValueError):stage_order(j)

    def test_finite_harmonic_growth_bound_and_zero_output_split(self):
        self.assertEqual(harmonic_bound(0),4);self.assertEqual(harmonic_bound(8),1024)
        self.assertEqual(max(abs(m) for m in next_harmonics([-4,-1,1,4])),8)
        self.assertNotIn(0,next_harmonics([-1,1]))

    def test_invalid_harmonic_sets_rejected(self):
        for v in ([],[0],[True],[F(1,2)]):
            with self.assertRaises(ValueError):next_harmonics(v)

    def test_physical_jet_loss_independent_of_stage(self):self.assertEqual(physical_loss(F(1,200),3),F(851,100))

    def test_each_fixed_physical_target_has_first_finite_stage(self):
        for m,target in ((0,0),(2,1),(5,10)):
            j=first_stage(target,m=m);self.assertGreaterEqual(physical_power(j,F(1,200),m),target)
            self.assertLess(physical_power(j-1,F(1,200),m),target)

    def test_diagnostic_nonnegative_zero_jet_needs_1002_unevaluated_stages(self):self.assertEqual(first_stage(0),1002)

    def test_low_stage_physical_residual_still_has_negative_power(self):self.assertLess(physical_power(0,F(1,200),0),0)

    def test_invalid_physical_domain_and_jet_order_rejected(self):
        for h,m in ((0,0),(F(1,100),0),(F(1,200),-1),(F(1,200),True)):
            with self.assertRaises(ValueError):physical_loss(h,m)

    def test_exact_all_twelve_complete_curl_divergence_components(self):
        r=harmonic_review();self.assertEqual(r['exact_complete_divergence_components'],12);self.assertTrue(r['checks']['all_complete_curl_harmonics_divergence_free'])

    def test_all_216_pairwise_longitudinal_reductions_including_zero_outputs(self):
        r=harmonic_review();self.assertEqual(r['exact_phase_reduction_components'],216)
        self.assertTrue(r['checks']['every_pair_phase_transport_reduces_without_half_frequency_loss']);self.assertTrue(r['checks']['zero_output_harmonics_do_not_require_division_by_output'])

    def test_complete_amplitudes_not_falsely_assumed_transverse(self):self.assertTrue(harmonic_review()['checks']['complete_amplitudes_not_assumed_exactly_transverse'])

    def test_higher_harmonic_source_reality_and_highest_mode_eight(self):
        c=harmonic_review()['checks'];self.assertTrue(c['conjugate_source_harmonics_preserve_reality']);self.assertTrue(c['highest_output_harmonic_eight_is_present'])

    def test_actual_mean_wave_crosses_and_new_squares_not_removed(self):
        c=harmonic_review()['checks'];self.assertTrue(c['mean_wave_transport_not_removed_by_wave_incompressibility']);self.assertTrue(c['new_wave_square_nonzero'])

    def test_all_90_full_nonlinear_increment_components(self):
        r=harmonic_review();self.assertEqual(r['exact_nonlinear_increment_components'],90);self.assertTrue(r['checks']['exact_full_nonlinear_increment_has_both_crosses_and_new_square'])

    def test_all_six_covariance_entries_and_independent_zero_modes(self):
        r=covariance_review();self.assertEqual(r['exact_zero_mode_covariance_components'],18);self.assertTrue(all(r['checks'].values()))

    def test_finite_diagnostics_not_actual_supported_global_fields(self):
        r=harmonic_review();self.assertTrue(r['graph_free_single_label_polynomial_diagnostic_not_actual_two_torus_pulse']);self.assertTrue(r['finite_endpoint_regularity_not_actual_C_infinity_certificate'])

    def test_all_943_historical_and_18_direct_source_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),943);self.assertEqual(len(self.record['source_hashes']),18)

    def test_frozen_record_and_JSON_roundtrip_match(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r);validate(json.loads(json.dumps(self.record)))

    def test_no_realization_stage_uniformity_or_blowup_promotion(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)

    def test_cached_reviews_not_mutable_aliases(self):
        r=generate();r['reviews']['finite_harmonic_complete_source_controls']['checks']['all_complete_curl_harmonics_divergence_free']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(harmonic_review()['checks']['all_complete_curl_harmonics_divergence_free'])

    def test_same_count_provenance_key_substitution_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])

    def test_protocol_gates_sources_remaining_not_erased(self):
        for key in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[key],dict):r[key].pop(next(iter(r[key])))
            else:r[key].pop()
            with self.assertRaises(ValueError):validate(r)

    def test_cli_refuses_overwriting_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')


if __name__=='__main__':unittest.main()
