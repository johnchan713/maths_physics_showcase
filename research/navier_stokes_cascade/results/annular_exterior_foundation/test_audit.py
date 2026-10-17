#!/usr/bin/env python3
"""Nonconstant annular maps, Gamma pole jets, reserved velocity/pressure and scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from math import factorial
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from exterior_bounds import heat_A,angular_majorant,derivative_bounds,foundation_review,tail_review,support_review
from exterior_checks import maps,annular_jets,triple_convolution_checks,heat_endpoint,heat_ode,reserved_plateau,gamma_stats,diagnostic_review
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,main


class ExteriorFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_parent_and_same_frequency(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'760ee0033020260296c1f6f3f3a264cfc0f4402c')
        self.assertEqual(self.record['protocol']['unchanged_frequency'],'N=1+floor(H^32)')

    def test_retained_core_and_finite_existence_dependencies_not_promoted(self):
        self.assertTrue(FLAGS['finite_order_existence_inverse_and_gap_lemmas_remain_explicit_dependencies'])
        self.assertFalse(FLAGS['all_order_weighted_T1_over_T0_and_background_edge_join_independently_verified'])

    def test_whole_compact_regular_input_gate_now_supplied_relatively(self):
        self.assertTrue(FLAGS['whole_compact_smooth_input_to_same_finite_N_maps_supplied_relative_to_exact_core_and_finite_existence_bounds'])
        self.assertTrue(all(foundation_review()['checks'].values()))

    def test_actual_annular_Q_eta_is_not_constant(self):
        self.assertIn('eta-dependent',foundation_review()['annular_Q'])
        self.assertEqual(foundation_review()['normalized_polynomial_derivative_suprema'],[4,8,8,4,1])

    def test_polynomial_normalized_derivative_norm_includes_every_nonzero_order(self):
        m=maps();eta=1+m.variable(0,n=1,order=8);p=(1+eta*eta)**2
        self.assertEqual([p.partial((k,))/factorial(k) for k in range(9)],[4,8,8,4,1,0,0,0,0])

    def test_nonzero_base_annular_root_matches_separate_jet_matrix_Newton(self):
        m=maps()
        for e in (Q(-1),Q(0),Q(1,3),Q(1)):
            c,n,z,eta,J=annular_jets(e)
            self.assertTrue(all(v.value()!=0 for v in c))
            self.assertTrue(all(m.zero(a-b) for a,b in zip(c,n)))
        self.assertEqual(diagnostic_review()['annular_coefficients'],180)

    def test_all_160_triple_derivative_components_keep_Q_eta_derivatives(self):
        for e in (Q(-1),Q(0),Q(1,3),Q(1)):
            c,n,z,eta,J=annular_jets(e);self.assertTrue(all(triple_convolution_checks(c,z,eta,J)))
        self.assertEqual(diagnostic_review()['triple_convolution_components'],160)

    def test_freezing_Q_eta_fails_even_at_symmetric_eta_zero(self):
        m=maps();c=annular_jets(Q(0))[0];wrong=annular_jets(Q(0),freeze_Q=True)[0]
        self.assertTrue(any(not m.zero(a-b) for a,b in zip(c,wrong)))

    def test_first_root_derivative_bound_keeps_Q_prime_times_c_squared(self):
        d=derivative_bounds(2,[0,3,4],[5,7,11],Q(1,10))
        self.assertEqual(d[1],4*(3+7*Q(1,10)**2))
        self.assertEqual(d[2],4*(4+5*d[1]**2+2*7*d[0]*d[1]+11*d[0]**2))

    def test_invalid_derivative_and_heat_parameters_rejected(self):
        for k in (-1,True,Q(1,2)):
            with self.assertRaises(ValueError):heat_A(Q(1,8),k)
        for h in (0,1,-1):
            with self.assertRaises(ValueError):heat_A(h,1)
        with self.assertRaises(ValueError):derivative_bounds(0,[1],[1],1)

    def test_exact_heat_Gamma_derivative_moments(self):
        h=Q(1,8);self.assertEqual(heat_A(h,0),1);self.assertEqual(heat_A(h,1),Q(9,64))
        for k in range(8):self.assertEqual(heat_A(h,k+1),heat_A(h,k)*(h+k)*(1+h+k))

    def test_all_18_heat_endpoint_chain_coefficients_at_both_poles(self):
        for e in (Q(-1),Q(1)):
            H,checks=heat_endpoint(e);self.assertTrue(all(checks))
        self.assertEqual(diagnostic_review()['heat_endpoint_coefficients'],18)

    def test_second_eta_heat_derivative_retains_Z_etaeta_term(self):
        h=Q(1,8);X=Q(64);H,_=heat_endpoint(Q(1),h,X)
        first_only=heat_A(h,2)*16/X**2
        self.assertEqual(H.partial((2,)),first_only+4*heat_A(h,1)/X)
        self.assertNotEqual(H.partial((2,)),first_only)

    def test_heat_differential_identity_all_available_exact_jets(self):
        self.assertTrue(maps().zero(heat_ode()))

    def test_Gamma_integral_matches_independent_Tricomi_U_for_ten_components(self):
        r=gamma_stats();self.assertEqual(r['count'],10);self.assertLess(r['maximum_error'],mp.mpf('1e-50'))
        self.assertTrue(r['bounds'])

    def test_all_fixed_angular_majorants_finite_at_poles_without_dividing_by_d(self):
        for k in (0,1,2,8,20):self.assertGreater(angular_majorant(Q(1,8),k),0)
        for e in (Q(-1),Q(1)):
            H,_=heat_endpoint(e)
            for k in range(1,9):self.assertLessEqual(abs(H.partial((k,))),angular_majorant(Q(1,8),k)/64)

    def test_heat_Taylor_coefficient_ratios_force_zero_analytic_radius(self):
        h=Q(1,8)
        for k in (1,2,8,1000):
            ratio=(k+h)*(k+1+h)/(k+1)
            self.assertGreater(ratio,k)
        self.assertFalse(FLAGS['heat_tail_analytic_in_eta_at_poles'])

    def test_positive_h_not_interchanged_with_infinite_integration(self):
        self.assertTrue(tail_review()['checks']['fixed_h_positive_kept_before_infinite_integration'])
        self.assertFalse(FLAGS['estimates_uniform_in_h_to_zero_or_all_orders_or_N'])

    def test_global_angular_moment_is_renormalized_not_raw_I(self):
        self.assertIn('renormalized_R=0',tail_review()['global_identities'])
        self.assertTrue(tail_review()['checks']['raw_I_total_not_set_to_zero_or_finite'])

    def test_exact_reserved_intervals_and_heat_compensation_gap(self):
        r=support_review();self.assertEqual(r['stage_offsets']['I3'],[-14,-9]);self.assertEqual(r['stage_offsets']['I4'],[-8,-3])
        self.assertTrue(all(r['checks'].values()))

    def test_reserved_leading_M_and_V_nonzero_despite_zero_axial_field(self):
        m=maps()
        for e in (Q(-1),Q(0),Q(1,3),Q(1)):
            a=reserved_plateau(e);self.assertTrue(m.zero(a['V']+3))
            self.assertEqual(a['M'].partial((0,1)),3)

    def test_reserved_pressure_uses_actual_density_and_full_K3_memory(self):
        m=maps();a=reserved_plateau()
        self.assertTrue(m.zero(a['Pi'].derivative(0)-a['E']**2/2))
        self.assertNotEqual(a['K'].value(),0);self.assertNotEqual(a['K'].partial((0,1)),0)

    def test_zero_M_shortcut_would_remove_required_radial_velocity(self):
        a=reserved_plateau(Q(0));self.assertEqual(a['M'].value(),0)
        self.assertEqual(a['V'].value(),-3) # value M=0 at one angle does not erase M_eta

    def test_leading_pressure_memory_not_equal_to_summed_mean_pressure(self):
        self.assertTrue(FLAGS['complete_intermediate_pressure_memory_and_exceptional_pressure_orders_retained'])
        self.assertIn('higher_pressure_tails_not_dropped',str(support_review()['checks']))

    def test_all_order_weighted_stress_quotient_is_separate_gate(self):
        self.assertTrue(support_review()['checks']['all_order_weighted_T1_over_T0_join_not_inferred_from_smoothness'])
        self.assertTrue(any('weighted' in x for x in self.record['remaining_obligations']))

    def test_uniform_physical_leading_bounds_not_promoted_to_actual_pulses(self):
        self.assertTrue(FLAGS['all_fixed_physical_leading_jet_bounds_on_finite_profile_ranges_joined'])
        self.assertFalse(FLAGS['unconditional_actual_pulse_bounds_verified'])

    def test_frozen_record_and_JSON_roundtrip_match_current_sources(self):
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.record,f)
        validate(json.loads(json.dumps(self.record)))

    def test_all_818_inherited_and_18_direct_bindings_preserved(self):
        self.assertEqual(len(self.record['inherited_provenance']),818)
        self.assertEqual(len(self.record['source_hashes']),18)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_cached_diagnostics_are_snapshotted_before_returning_records(self):
        r=generate();key='freezing_Q_eta_misses_derivatives'
        r['reviews']['independent_nonconstant_root_and_heat_controls']['checks'][key]=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(diagnostic_review()['checks'][key])

    def test_inherited_provenance_cannot_be_mutated_or_substituted(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['invented:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('invented:key',parent_review()[1])

    def test_unconditional_or_analyticity_scope_promotions_rejected(self):
        for key,value in FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):validate(r)

    def test_pressure_or_reserved_plateau_formula_tampering_rejected(self):
        for key in ('pressure','plateau'):
            r=deepcopy(self.record);r['reviews']['actual_reserved_plateau_and_pressure_memory'][key]='0'
            with self.assertRaises(ValueError):validate(r)

    def test_remaining_gate_or_binding_cannot_be_erased(self):
        for key in ('remaining_obligations','gates','source_hashes'):
            r=deepcopy(self.record)
            if key=='remaining_obligations':r[key].pop()
            else:r[key].pop(next(iter(r[key])))
            with self.assertRaises(ValueError):validate(r)

    def test_cli_refuses_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':unittest.main()
