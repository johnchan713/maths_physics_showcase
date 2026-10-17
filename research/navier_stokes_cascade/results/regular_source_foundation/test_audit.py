#!/usr/bin/env python3
"""Regular linear/source bounds with exact-jet and explicit scope controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from source_bounds import KAPPA,linear_budget,nonlinear_budget,initial_budget,next_harmonics,review
from source_checks import linear_checks,transport_checks,endpoint_checks
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main


class RegularSourceFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_live_parent_and_frozen_parameters(self):
        self.assertEqual(protocol()['parent_commit'],'10b5bdff358cd6c54474be533dec1618f829d940')
        self.assertIn('t1=C^-20',protocol()['unchanged_parameters'])

    def test_sharper_linear_gain_exact(self):
        for a in (F(-7),F(0),F(1,2),F(8)):
            r=linear_budget(a);self.assertEqual(r['target']-a,F(1,2)-KAPPA)
            self.assertTrue(all(x>=r['target'] for x in r['terms'].values()))

    def test_old_three_kappa_gain_still_valid(self):
        self.assertGreater(F(1,2)-KAPPA,F(1,2)-3*KAPPA)
        self.assertTrue(review()['checks']['former_three_kappa_bound_remains_valid'])

    def test_curl_fast_time_does_not_charge_another_graph_loss(self):
        r=linear_budget(2);self.assertEqual(r['terms']['curl_fast_time'],r['curl_remainder'])
        self.assertEqual(r['terms']['curl_K'],r['target'])
        self.assertEqual(r['terms']['curl_principal_damping'],r['target'])

    def test_pressure_and_curl_remainder_orders_not_conflated(self):
        r=linear_budget(2);self.assertEqual(r['pressure']-r['curl_remainder'],KAPPA)

    def test_remaining_curl_viscosity_is_stronger(self):
        r=linear_budget(2);self.assertGreater(r['terms']['curl_remaining_viscosity'],r['target'])

    def test_both_pressure_and_mixed_viscosity_limit_gain(self):
        r=linear_budget(1)
        for k in ('pressure_gradient','mixed_viscosity_and_phase_divergence'):self.assertEqual(r['terms'][k],r['target'])

    def test_invalid_kappa_rejected(self):
        for k in (0,F(-1),F(1,4),1):
            with self.assertRaises(ValueError):linear_budget(1,k)
            with self.assertRaises(ValueError):nonlinear_budget(1,1,1,k)

    def test_all_672_exact_full_linear_jet_components(self):
        r=linear_checks();self.assertEqual(r['exact_linear_jet_components'],672);self.assertTrue(r['principal_plus_remainder_matches'])

    def test_angular_viscous_connection_cannot_be_omitted(self):self.assertTrue(linear_checks()['omission_controls']['angular'])
    def test_mixed_phase_amplitude_viscosity_cannot_be_omitted(self):self.assertTrue(linear_checks()['omission_controls']['mixed'])
    def test_pressure_amplitude_gradient_cannot_be_omitted(self):self.assertTrue(linear_checks()['omission_controls']['pressure'])
    def test_repeated_radial_graph_differentiates_its_multiplier(self):self.assertTrue(linear_checks()['omission_controls']['graph_coefficient'])

    def test_all_210_exact_complete_transport_jet_components(self):
        r=transport_checks();self.assertEqual(r['exact_transport_jet_components'],210);self.assertTrue(r['complete_harmonic_transport_matches'])

    def test_complete_curl_longitudinal_identity(self):self.assertTrue(transport_checks()['exact_complete_curl_longitudinal_identity'])
    def test_complete_amplitude_is_not_principal_transverse_amplitude(self):self.assertTrue(transport_checks()['complete_amplitude_not_transverse'])

    def test_zero_and_negative_output_harmonics_are_retained(self):
        r=transport_checks();self.assertIn(0,r['output_multipliers']);self.assertIn(-3,r['output_multipliers'])

    def test_same_label_cancellation_recovers_half_power_minus_kappa(self):
        r=nonlinear_budget(F(1,2),F(3,2),1)
        self.assertEqual(r['same_label_wave_wave']-r['apparent_uncancelled_wave_phase'],F(1,2)-KAPPA)

    def test_mean_wave_half_power_loss_is_not_a_graph_loss(self):
        r=nonlinear_budget(F(1,2),F(1,2),1)
        self.assertEqual(r['mean_wave'],1)
        self.assertLess(r['tangential_mean_phase'],r['wave_on_mean'])

    def test_different_phase_and_mean_wave_frequency_witnesses(self):
        # Differentiate and convolve sparse Fourier fields exactly, then use
        # Parseval. Coefficients are Gaussian rationals, not floating samples.
        def mul(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
        def advect(a,b):
            out={}
            for u,aa in a.items():
                for v,bb in b.items():
                    mode=tuple(x+y for x,y in zip(u,v));out.setdefault(mode,[[F(0),F(0)] for _ in range(3)])
                    contraction=(-sum(x[1]*j for x,j in zip(aa,v)),sum(x[0]*j for x,j in zip(aa,v)))
                    for i,z in enumerate(bb):
                        x=mul(contraction,z);out[mode][i]=[p+q for p,q in zip(out[mode][i],x)]
            return out
        def field(k,axis,component):
            out={}
            for sign in (-1,1):
                mode=tuple(sign*k if i==axis else 0 for i in range(3))
                out[mode]=[(F(1,2) if i==component else F(0),F(0)) for i in range(3)]
            return out
        norm=lambda d:sum(re*re+im*im for row in d.values() for re,im in row)
        for k in (2,8,32):
            a,b=field(k,0,1),field(k,1,0);ab,ba=advect(a,b),advect(b,a)
            combined={mode:[[x+y for x,y in zip(v,w)] for v,w in zip(ab[mode],ba[mode])] for mode in ab}
            self.assertEqual(norm(combined),F(k*k,2))
            mean={(0,0,0):[(F(0),F(0)),(F(0),F(0)),(F(1),F(0))]}
            self.assertEqual(norm(advect(mean,field(k,2,0))),F(k*k,2))

    def test_actual_initial_source_all_three_groups_one_minus_kappa(self):
        self.assertEqual(set(initial_budget().values()),{1-KAPPA})

    def test_harmonic_grammar_excludes_zero_from_forward_sources(self):
        current={-2,-1,1,2};out=next_harmonics(current)
        self.assertNotIn(0,out);self.assertEqual(max(map(abs,out)),4)
        for _ in range(5):
            out=next_harmonics(out);self.assertEqual(min(out),-max(out))

    def test_invalid_forward_harmonics_rejected(self):
        for v in ([],[0],[True],[F(1)]):
            with self.assertRaises(ValueError):next_harmonics(v)

    def test_14_endpoint_jet_controls_are_not_C_infinity_certificate(self):
        r=endpoint_checks();self.assertEqual(r['endpoint_jet_controls'],14);self.assertTrue(r['cut_potential_and_pressure_endpoint_jets'])
        self.assertIn('C7 only',r['fixture_regularity'])

    def test_flat_tail_witness_is_not_forward_inverted(self):
        r=endpoint_checks();self.assertIn('core response nonflat',r['flat_homogeneous_tail_witness'])
        self.assertTrue(r['flat_tails_must_remain_additive_not_forward_sources'])

    def test_whole_pressure_and_background_tails_not_erased(self):self.assertTrue(FLAGS['complete_Pi1_Pi2_background_pressure_and_mean_moments_retained'])

    def test_full_iteration_and_sharper_cycle_not_promoted(self):
        for key,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):validate(r)

    def test_all_889_historical_and_18_direct_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),889);self.assertEqual(len(self.record['source_hashes']),18)

    def test_frozen_record_and_JSON_roundtrip_match(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r)
        validate(json.loads(json.dumps(self.record)))

    def test_cached_exact_reviews_are_not_mutable_record_aliases(self):
        r=generate();r['reviews']['exact_full_viscous_linear_expansion']['principal_plus_remainder_matches']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(linear_checks()['principal_plus_remainder_matches'])

    def test_provenance_key_substitution_rejected_even_at_same_count(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])

    def test_formula_gates_sources_and_remaining_obligations_cannot_be_erased(self):
        for key in ('protocol','gates','source_hashes','remaining_obligations'):
            r=deepcopy(self.record)
            if isinstance(r[key],dict):r[key].pop(next(iter(r[key])))
            else:r[key].pop()
            with self.assertRaises(ValueError):validate(r)

    def test_cli_refuses_overwriting_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'evidence.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')


if __name__=='__main__':unittest.main()
