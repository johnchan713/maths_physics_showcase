#!/usr/bin/env python3
"""Closed-time geometry, exact heat/pressure jets, endpoint limits and scope."""
from copy import deepcopy
from fractions import Fraction as Q
from math import factorial
import json
from pathlib import Path
import tempfile
import unittest
from endpoint import (floor_log2,ceil_log2,band_window,active_stage_bound,closed_centers,
    requirements,heat_A,heat_S,normalized_heat_time_jet,normalized_pressure_time_jet,
    endpoint_review,nonanalytic_review)
from endpoint_checks import q_value,chart_stats,heat_stats,Cauchy_stats,diagnostic_review
from audit import HERE,FLAGS,generate,validate,compare_records,main,parent_review


class CurrentEndpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_input_is_current_realized_field(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'aac13c80650aa0f46320e723270431b9d2d64b3f')
        self.assertTrue(self.record['gates']['current_actual_realization_replayed'])

    def test_historical_force_conclusion_not_imported_for_current_field(self):
        r=self.record['reviews']['historical_endpoint_transfer_replayed']
        self.assertTrue(r['historical_force_and_energy_conclusions_not_imported_for_current_field'])
        self.assertTrue(all(r['checks'].values()))

    def test_closed_band_window_includes_boundary_bands(self):
        self.assertEqual(band_window(Q(2)**-24,Q(2)**-20),dict(first=19,last=25,count=7))

    def test_nondyadic_band_window_is_exact_not_float_log(self):
        self.assertEqual(band_window(Q(3,2**25),Q(3,2**21)),dict(first=19,last=24,count=6))

    def test_fixed_first_band_can_leave_empty_window(self):
        self.assertEqual(band_window(Q(2)**-24,Q(2)**-20,100)['count'],0)

    def test_arbitrarily_small_positive_dyadic_logs_are_exact(self):
        x=Q(2)**-100000
        self.assertEqual(floor_log2(x),-100000);self.assertEqual(ceil_log2(x),-100000)
        self.assertEqual(floor_log2(3*x),-99999);self.assertEqual(ceil_log2(3*x),-99998)

    def test_invalid_band_domains_and_logs_rejected(self):
        for args in ((0,Q(1,2)),(Q(1,2),Q(1,4)),(Q(1,2),1),(Q(1,4),Q(1,2),True)):
            with self.assertRaises(ValueError):band_window(*args)
        for f in (floor_log2,ceil_log2):
            for x in (0,-1):
                with self.assertRaises(ValueError):f(x)

    def test_closed_terminal_time_mesh_indices_are_already_selected(self):
        for ell in (16,100,10**6):self.assertEqual(closed_centers(0,0,Q(1,ell**6)),(-1,1))

    def test_negative_and_nondyadic_mesh_edges_use_ceil_floor_correctly(self):
        self.assertEqual(closed_centers(Q(-3,2),Q(3,2),1),(-2,2))
        self.assertEqual(closed_centers(-1,1,1),(-2,2))

    def test_invalid_mesh_intervals_rejected(self):
        for args in ((1,0,1),(0,1,0),(0,1,-1)):
            with self.assertRaises(ValueError):closed_centers(*args)

    def test_cutoff_stage_bound_keeps_strict_power_of_two_inequality(self):
        self.assertEqual(active_stage_bound(Q(2)**-24,Q(2)**-16),7)
        self.assertEqual(active_stage_bound(Q(1,4),Q(1,2)),0)

    def test_nondyadic_stage_bound_does_not_drop_last_possible_stage(self):
        self.assertEqual(active_stage_bound(Q(3,2**26),Q(2)**-16),8)

    def test_invalid_stage_bound_domain_rejected(self):
        for args in ((0,Q(1,2)),(Q(1,2),Q(1,4)),(Q(1,2),1)):
            with self.assertRaises(ValueError):active_stage_bound(*args)

    def test_all_fixed_terminal_and_future_force_derivative_requirements(self):
        for m in (0,2,10000):
            r=requirements(m)
            self.assertEqual(r['representative_Cauchy_bound'],m+1)
            self.assertEqual(r['velocity_Cauchy_potential_bound'],m+2)
            self.assertEqual(r['future_force_Cauchy_potential_bound'],m+4)
            self.assertEqual(r['future_force_Cauchy_pressure_bound'],m+2)
        for m in (-1,True,1.5):
            with self.assertRaises(ValueError):requirements(m)

    def test_naive_half_domain_does_not_ensure_enlarged_chart_margin(self):
        q,Qchart,qstar=Q(1,4),Q(1,2),Q(1)
        self.assertLess(q,qstar/2);self.assertEqual(q/Qchart,Q(1,2));self.assertGreater(4*Qchart,qstar)
        self.assertTrue(endpoint_review()['checks']['all_active_enlarged_charts_stay_inside_background_domain'])

    def test_punctured_cover_has_strict_heat_exterior_margin(self):
        self.assertGreater(Q(16,3),1)
        self.assertTrue(endpoint_review()['checks']['uniform_axis_hole_for_annular_fields_from_Xa_and_positive_q'])

    def test_terminal_q_at_both_poles_and_singular_corner_rejected(self):
        for z in (-.3,.3):self.assertAlmostEqual(q_value(z,1),abs(z)**(1/.495))
        for z,t in ((0,1),(.3,1.01)):
            with self.assertRaises(ValueError):q_value(z,t)

    def test_terminal_implicit_q_all_third_order_jets_and_time_sign(self):
        r=chart_stats();self.assertLess(r['implicit_equation_error'],2e-12)
        self.assertLess(r['q_time_sign_error'],2e-12)
        self.assertLess(r['slow_time_chart_error'],2e-12)
        self.assertGreater(r['wrong_slow_time_power_gap'],1e-3)

    def test_complete_tuple_and_velocity_jet_chart_gluing(self):
        r=chart_stats();self.assertEqual(r['chart_comparisons'],12);self.assertLess(r['chart_jet_error'],2e-10)

    def test_wrong_independent_fixed_R_bumps_do_not_glue(self):
        self.assertGreater(chart_stats()['wrong_fixed_chart_bump_gap'],1e-3)

    def test_terminal_curl_and_direct_swirl_remain_solenoidal(self):
        self.assertLess(chart_stats()['divergence_error'],2e-12)

    def test_zero_extending_nonzero_terminal_pressure_would_jump(self):
        self.assertGreater(chart_stats()['zero_time_extension_pressure_jump'],1e-3)

    def test_independent_Gamma_integrals_match_heat_time_jets(self):
        self.assertLess(heat_stats()['heat_moment_error'],1e-45)

    def test_full_pressure_convolution_matches_independent_moments(self):
        self.assertLess(heat_stats()['complete_pressure_convolution_error'],1e-45)

    def test_exact_first_heat_and_pressure_jets_retain_time_sign_and_memory(self):
        h=Q(1,200)
        self.assertEqual(normalized_heat_time_jet(h,0),1)
        self.assertEqual(normalized_heat_time_jet(h,1),Q(201,10000))
        self.assertEqual(normalized_pressure_time_jet(h,0),Q(-50,101))
        self.assertEqual(normalized_pressure_time_jet(h,1),Q(-1,100))

    def test_pressure_second_jet_cannot_drop_square_cross(self):
        h=Q(1,200)
        self.assertEqual(heat_S(h,2)-2*heat_A(h,2),2*heat_A(h,1)**2)
        self.assertGreater(heat_S(h,2),2*heat_A(h,2))

    def test_heat_Taylor_ratio_identity_grows_without_bound(self):
        self.assertTrue(all(nonanalytic_review()['checks'].values()))
        for j in (1,1000):
            h=Q(1,200);ratio=normalized_heat_time_jet(h,j+1)/normalized_heat_time_jet(h,j)/(j+1)
            self.assertEqual(ratio-4*j,4*h*(2*j+1+h)/(j+1))

    def test_factorial_obstruction_survives_arbitrarily_small_fixed_positive_h(self):
        for h in (Q(1,200),Q(1,10**30)):
            for j in (1,2,40):
                self.assertGreaterEqual(heat_A(h,j),h*factorial(j-1)*factorial(j))
                self.assertGreaterEqual(heat_S(h,j),2*heat_A(h,j))

    def test_h_zero_or_invalid_jet_order_cannot_erase_nonanalytic_obstruction(self):
        for h,j in ((0,1),(Q(1,100),1),(Q(1,200),-1),(Q(1,200),True)):
            with self.assertRaises(ValueError):normalized_heat_time_jet(h,j)

    def test_finite_Taylor_remainder_has_uniform_known_bound(self):
        self.assertLess(Cauchy_stats()['maximum_finite_Taylor_remainder_ratio'],1)

    def test_wrong_tau_time_sign_fails_finite_Taylor_estimate(self):
        self.assertGreater(Cauchy_stats()['reversed_time_Taylor_gap'],100)

    def test_bounded_values_and_origin_flatness_are_not_enough(self):
        self.assertGreater(min(Cauchy_stats()['bounded_value_without_next_derivative_limit_gaps']),1.99)

    def test_exact_Gaussian_exclusion_and_complete_current_seed_preserved(self):
        self.assertTrue(parent_review()[0]['flags']['designated_Gaussian_tails_never_inverted_or_independently_summed'])
        self.assertTrue(FLAGS['actual_flat_velocities_base_pressure_tails_and_Gaussian_exclusion_preserved'])

    def test_716_inherited_hashes_and_16_source_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),716)
        self.assertTrue(all(self.record['inherited_provenance'].values()))
        self.assertEqual(len(self.record['source_hashes']),16)

    def test_frozen_evidence_reproduces(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r)

    def test_record_roundtrips_without_numpy_or_near_zero_noise_refreezing(self):
        r=json.loads(json.dumps(self.record));validate(r);compare_records(self.record,r)
        d=self.record['reviews']['manufactured_terminal_and_failure_controls']
        self.assertTrue(d['raw_roundoff_scores_not_refrozen_under_an_ancestor_comparator'])

    def test_unconditional_force_origin_regular_and_analytic_promotions_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError):validate(r)

    def test_inherited_uniform_physical_premises_cannot_be_erased(self):
        r=deepcopy(self.record)
        r['flags']['inherited_all_order_physical_background_and_inverse_bounds_remain_analytical_premises']=False
        with self.assertRaises(ValueError):validate(r)

    def test_changed_pressure_jet_or_missing_controls_rejected_despite_passed_gates(self):
        for name in ('pressure','control'):
            r=deepcopy(self.record)
            if name=='pressure':r['reviews']['exact_exterior_swirl_and_pressure_jets']['rows'][1]['normalized_pressure_jet']='0'
            else:r['reviews']['manufactured_terminal_and_failure_controls']['checks'].pop('independent_fixed_R_pressure_bumps_do_not_glue')
            with self.assertRaises(ValueError):validate(r)

    def test_missing_gate_obligation_or_provenance_rejected(self):
        for section in ('gates','remaining_obligations','source_hashes','inherited_provenance'):
            r=deepcopy(self.record)
            if section=='remaining_obligations':r[section].pop()
            else:r[section].pop(next(iter(r[section])))
            with self.assertRaises(ValueError):validate(r)

    def test_cli_does_not_overwrite_frozen_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':unittest.main()
