#!/usr/bin/env python3
"""All-scale primary physical-pulse identities and explicit remaining gates."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from pulse_bounds import exponents,harmonic_loss,primary_powers,riccati_outward,tail_square_constant,input_review,pulse_review,physical_review,support_review
from pulse_checks import maps,band_jets,graph_stats,frame_stats,scalar_stats,diagnostic_review
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,main


class PhysicalPulseFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_actual_background_parent_and_frozen_radius(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'585342bc4225cd69cf07a0a9a24d65ef01729804')
        self.assertIn('r0=1/6057072',self.record['protocol']['unchanged_parameters'])

    def test_exact_core_and_finite_lemma_dependencies_not_erased(self):
        self.assertTrue(FLAGS['exact_core_and_finite_existence_inverse_gap_lemmas_remain_dependencies'])
        self.assertFalse(FLAGS['unconditional_all_actual_pulse_families_verified'])

    def test_enlarged_chart_controls_q_not_only_Q(self):
        self.assertIn('4Q<qstar',self.record['protocol']['enlarged_chart'])
        self.assertLess(Q(3,4),1);self.assertGreater(4*Q(3,4),1)

    def test_implicit_band_jets_including_one_sided_terminal_points(self):
        h=Q(1,200)
        for e in (Q(-1),Q(1,3),Q(1)):
            s,r=band_jets(e);self.assertFalse(r.terms)
            L=1-2*h*e*e
            self.assertEqual(s.partial((0,1)),1/L)
            self.assertEqual(s.partial((1,0)),2*e/L)
        self.assertEqual(diagnostic_review()['implicit_band_components'],84)

    def test_no_Q_power_in_implicit_coordinate_equation(self):
        for h in (Q(1,200),Q(9,1000)):
            e=exponents(h);self.assertEqual(2*e['D']+2*h,1)

    def test_representatives_and_discrete_carriers_frozen_before_derivatives(self):
        self.assertTrue(self.record['protocol']['discrete_data_frozen_before_all_differentiation'])
        self.assertTrue(input_review()['checks']['representatives_selected_inside_closed_active_chart_intersection'])

    def test_nearest_nonzero_rounding_may_exceed_half_over_k(self):
        self.assertGreater(abs(1-Q(1,10)),Q(1,2))
        self.assertLessEqual(abs(1-Q(1,10)),1)
        self.assertTrue(input_review()['checks']['entire_pulse_L_over_k_rounding_and_fixed_radius_costs_retained'])

    def test_all_420_independent_mixed_physical_curl_pressure_components(self):
        r=graph_stats();self.assertEqual(r['count'],420);self.assertTrue(r['matches'])

    def test_direct_complete_cartesian_curls_are_divergence_free(self):
        self.assertTrue(graph_stats()['divergences'])

    def test_radial_auxiliary_graph_cancels_for_all_checked_phase_jets(self):
        self.assertTrue(graph_stats()['phase_matches'])
        self.assertTrue(graph_stats()['radial_graph_omission_nonzero'])

    def test_amplitude_radial_graph_and_fast_time_cannot_be_omitted(self):
        r=graph_stats();self.assertTrue(r['radial_graph_omission_nonzero']);self.assertTrue(r['fast_time_omission_nonzero'])

    def test_repeated_graph_derivatives_must_differentiate_the_radial_coefficient(self):
        self.assertTrue(graph_stats()['frozen_graph_coefficient_omission_nonzero'])

    def test_cylindrical_connection_is_nonzero_in_complete_curl(self):
        self.assertTrue(graph_stats()['cylindrical_connection_nonzero'])

    def test_all_300_reference_frame_projection_left_inverse_components(self):
        r=frame_stats();self.assertEqual(r['reference_components'],300);self.assertTrue(r['reference'])

    def test_all_216_moving_frame_jets_retain_inverse_length_scaling(self):
        r=frame_stats();self.assertEqual(r['moving_components'],216);self.assertTrue(r['moving_scaling'])

    def test_pressure_preserves_constraint_with_moving_and_source_normal_terms(self):
        r=frame_stats();self.assertTrue(r['normal_pressure']);self.assertTrue(r['omitted_moving_or_source_normal_detected'])
        self.assertEqual(r['pressure_components'],12)

    def test_normal_floor_is_not_an_angular_component_floor(self):
        self.assertTrue(pulse_review()['checks']['tangential_normal_floor_not_inverse_of_angular_component_used'])
        self.assertIn('nprime.t',pulse_review()['pressure'])

    def test_all_350_phase_partition_costs_and_primary_normalizations(self):
        r=scalar_stats();self.assertEqual(r['partitions'],350);self.assertTrue(r['dominance']);self.assertTrue(r['normalizations'])

    def test_exact_primary_velocity_pressure_potential_physical_powers(self):
        h=Q(1,200);e=exponents(h);p=primary_powers(h,0,0)
        self.assertEqual(p,dict(velocity=-e['sx'],pressure=-2*e['sx'],potential=0))
        loss,upper=harmonic_loss(h,2,3,5);self.assertEqual(loss,upper)

    def test_invalid_h_orders_or_empty_phase_partitions_rejected(self):
        for h in (0,Q(1,100),-1):
            with self.assertRaises(ValueError):exponents(h)
        for a,b,j in ((0,0,0),(1,0,2),(-1,0,1),(True,0,1)):
            with self.assertRaises(ValueError):harmonic_loss(Q(1,200),a,b,j)

    def test_27_Riccati_barriers_have_strict_outward_sign(self):
        r=scalar_stats();self.assertEqual(r['barrier_components'],27);self.assertTrue(r['barriers'])
        ratio,out,upper=riccati_outward(1,1,8)
        self.assertEqual(ratio,Q(1,8));self.assertEqual(out,upper);self.assertLess(out,0)

    def test_short_pulse_rejected_before_lower_bound_claim(self):
        with self.assertRaises(ValueError):riccati_outward(1,1,7)
        self.assertFalse(FLAGS['positive_pointwise_lower_bound_for_oscillating_cutoff_velocity'])

    def test_actual_Haar_area_and_fixed_geometry_radius_costs_retained(self):
        self.assertIn('r0^(3/2)',pulse_review()['radius_weights'])
        self.assertIn('r0^(-3/4)',pulse_review()['radius_weights'])

    def test_all_243_exact_Gaussian_square_completions_keep_same_domain(self):
        r=scalar_stats();self.assertEqual(r['square_components'],243);self.assertTrue(r['squares'])
        d,c=tail_square_constant(Q(1,1000),20,20,100)
        self.assertEqual(c,d*d/(4*Q(1,1000)))

    def test_tail_constant_is_not_a_derivative_dependent_start_band(self):
        self.assertTrue(FLAGS['same_domain_Gaussian_tail_flatness_constants_derived_for_every_fixed_order'])
        self.assertIn('SAME entire band domain',support_review()['tail'])
        with self.assertRaises(ValueError):tail_square_constant(0,1,1,1)
        with self.assertRaises(ValueError):tail_square_constant(1,1,1,-1)

    def test_exact_square_root_edge_constants_and_closed_supports_retained(self):
        r=support_review();self.assertEqual(r['square_root_weight'],'exp(-t1^2/(2*ya^2)-2/yb^2)')
        self.assertEqual(r['supports']['pulse'],'abs(v/Ls-1/2)<=2/5')

    def test_positive_P_endpoints_do_not_supply_source_zero_extension(self):
        r=support_review()['checks']
        self.assertTrue(r['positive_endpoint_envelope_P_does_not_supply_zero_extension'])
        self.assertTrue(r['source_requires_its_own_full_harmonic_endpoint_extension'])

    def test_actual_primary_bounds_not_promoted_to_all_generated_sources(self):
        self.assertTrue(FLAGS['uniform_primary_complete_curl_pressure_and_fixed_potential_bounds_derived'])
        self.assertFalse(FLAGS['generated_nonlinear_source_bounds_and_supports_all_verified'])
        self.assertTrue(self.record['protocol']['source_mapping_not_generated_source_verification'])

    def test_complete_means_keep_radial_axial_velocity_and_pressure_tails_open(self):
        self.assertTrue(FLAGS['complete_background_radial_flow_and_I4_pressure_tails_retained'])
        self.assertFalse(FLAGS['unconditional_complete_covariance_mean_corrections_verified'])
        self.assertTrue(any('Pi1_Pi2' in x for x in self.record['remaining_obligations']))

    def test_frozen_record_and_JSON_roundtrip_match_sources(self):
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.record,f)
        validate(json.loads(json.dumps(self.record)))

    def test_all_855_inherited_and_19_direct_source_bindings_preserved(self):
        self.assertEqual(len(self.record['inherited_provenance']),855)
        self.assertEqual(len(self.record['source_hashes']),19)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_cached_reviews_are_snapshots_not_mutable_canonical_evidence(self):
        r=generate();key='complete_normal_pressure_retains_moving_and_source_terms'
        r['reviews']['independent_exact_physical_frame_and_tail_controls']['checks'][key]=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(diagnostic_review()['checks'][key])

    def test_provenance_key_substitution_rejected_even_with_same_count(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])

    def test_unconditional_positive_velocity_or_numerical_scope_promotions_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)

    def test_physical_pressure_or_endpoint_formula_tampering_rejected(self):
        for group,key in (('moving_frame_positive_pulse_and_Haar_weights','pressure'),('exact_supports_and_same_domain_Gaussian_tails','tail')):
            r=deepcopy(self.record);r['reviews'][group][key]='0'
            with self.assertRaises(ValueError):validate(r)

    def test_remaining_gates_bindings_and_source_obligations_not_erased(self):
        for key in ('remaining_obligations','gates','source_hashes'):
            r=deepcopy(self.record)
            if key=='remaining_obligations':r[key].pop()
            else:r[key].pop(next(iter(r[key])))
            with self.assertRaises(ValueError):validate(r)

    def test_cli_refuses_overwriting_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')


if __name__=='__main__':unittest.main()
