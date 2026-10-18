"""Canonical seed, finite operators, current-target and quantifier regressions."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,protocol,generate,validate,compare_records,main,initial_record,complete_controls
from cycle_bounds import K,operator_budget,prefix_orders,prefix_review,input_review,signed_coefficients,signed_review,cycle_review,physical_review,diagonal_review
from freshness import State,Input,build,review as freshness_review


class CanonicalCycleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()
    def test_current_parent_and_protocol(self):
        self.assertEqual(protocol()['parent_commit'],'5e49006bab195749c75635a710e1d13216ee7f2b')
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),protocol())
    def test_current_canonical_seed_inputs(self):self.assertTrue(all(input_review(initial_record())['checks'].values()))
    def test_wrong_initial_order_cannot_supply_new_seed(self):
        r=deepcopy(initial_record());r['reviews']['canonical_initial_order_ledger']['orders']['final_defects']='1'
        self.assertFalse(all(input_review(r)['checks'].values()))
    def test_fresh_prefix_computed_orders(self):
        p=prefix_orders();self.assertEqual(p['final_regular'],F(3,2)-2*K)
        self.assertEqual(p['final_full_mean'],p['final_regular']+F(1,2))
        self.assertTrue(all(prefix_review()['checks'].values()))
    def test_generic_covariance_not_real_primary_shortcut(self):
        p=prefix_orders();self.assertEqual(p['generic_signed_Haar'],2-2*K)
        self.assertEqual(p['real_signed_averaged'],3-6*K)
        self.assertFalse(FLAGS['initial_real_covariance_shortcut_used_for_complex_cycle'])
    def test_old_error_stays_at_rho(self):
        p=prefix_orders();self.assertEqual(p['rho'],1-K);self.assertLess(p['rho'],p['final_regular'])
    def test_exact_arbitrary_signed_inverse(self):self.assertTrue(all(signed_review()['checks'].values()))
    def test_signed_target_does_not_need_positive_updated_amplitude(self):
        a=[F(1,2),F(2,3)];da=signed_coefficients([[1,1],[-1,1]],a,[-100,300],F(1,200))
        self.assertTrue(any(x+y<0 for x,y in zip(a,da)))
    def test_zero_fixed_root_and_zero_epsilon_rejected(self):
        for a,eps in (([0,1],1),([1,1],0),([-1,1],1)):
            with self.assertRaises(ValueError):signed_coefficients([[1,1],[-1,1]],a,[1,2],eps)
    def test_singular_H_rejected(self):
        with self.assertRaises(ValueError):signed_coefficients([[1,1],[1,1]],[1,1],[1,2],1)
    def test_whole_half_line_affine_budget(self):self.assertTrue(all(cycle_review()['checks'].values()))
    def test_larger_available_gain_is_blocked(self):
        with self.assertRaises(ValueError):operator_budget().certificate(F(1,2))
    def test_negative_affine_slope_cannot_pass_by_boundary(self):
        o=operator_budget()
        with self.assertRaises(ValueError):o.line(0,100).margin(o.line(1,0))
    def test_fast_viscosity_and_stale_target_gaps(self):
        c=operator_budget().certificate()
        self.assertEqual(c['margins']['fast_full_mean:radial_viscosity'],F(1,2)-2*K)
        self.assertEqual(1-(operator_budget().MU-K),2*K)
    def test_stage_orders_are_finite_and_increase(self):
        o=operator_budget();self.assertEqual(o.stage_order(4)-o.stage_order(3),F(1,2)-K)
        for j in (-1,True,F(1,2)):
            with self.assertRaises(ValueError):o.stage_order(j)
    def test_every_fixed_physical_target_has_finite_stage(self):self.assertTrue(all(physical_review()['checks'].values()))
    def test_seed_physical_power_is_not_positive(self):self.assertLess(operator_budget().physical_power(0,F(1,200),0),0)
    def test_finite_harmonic_rule_keeps_zero_output_separate(self):
        o=operator_budget();self.assertNotIn(0,o.next_harmonics([-2,-1,1,2]))
        self.assertEqual(o.harmonic_bound(3),32)
    def test_all_six_complete_fourier_covariances(self):
        self.assertTrue(all(complete_controls().covariance_review()['checks'].values()))
    def test_full_nonlinear_increment_and_mean_transport(self):
        r=complete_controls().harmonic_review();self.assertTrue(all(r['checks'].values()))
        self.assertIn(0,r['nonzero_source_harmonics'])
    def test_current_bundle_prefix_and_cycles(self):
        for j in (0,1,4):
            s=build(j);self.assertTrue(s.clean());self.assertEqual(s.bundle,'canonical-reference-roots')
            self.assertEqual(s.wave,3+2*j);self.assertEqual(s.mean,6+2*j)
    def test_parent_tail_ledger_is_retained(self):
        s=build(2);self.assertEqual(s.tails[0],('parent-existing-ledger',0));self.assertEqual(len(s.tails),8)
    def test_missing_pressure_recompute_rejected(self):
        s=State();slow=s.read('slow');s.apply(s.read('fast'),'fast')
        with self.assertRaises(ValueError):s.apply(slow,'slow')
    def test_stale_prefast_target_rejected_even_after_recompute(self):
        s=State();slow=s.read('slow');s.apply(s.read('fast'),'fast');s.recompute()
        with self.assertRaises(ValueError):s.apply(slow,'slow')
    def test_historical_root_bundle_mixing_rejected(self):
        with self.assertRaises(ValueError):State().apply(Input('historical-roots',(0,2),'signed'),'signed')
    def test_additive_gaussian_tail_cannot_be_inverse_source(self):
        s=State()
        with self.assertRaises(ValueError):s.apply(s.read('tail'),'sourced')
    def test_source_stamps_are_only_a_protocol_model(self):self.assertTrue(all(freshness_review()['checks'].values()))
    def test_arbitrary_diagonal_counterexample(self):
        r=diagonal_review();self.assertTrue(all(r['checks'].values()))
        for row in r['rows']:
            self.assertGreater(F(row['fixed_stage_power']),0);self.assertGreater(F(row['log2_diagonal_residual']),0)
    def test_full_pressure_formula_and_moving_memory_remain(self):
        self.assertIn('nprime.t_m+n.f_m',protocol()['inverse'])
        self.assertIn('Pi1/Pi2',protocol()['pressure']);self.assertIn('whole partial_Z',protocol()['pressure'])
    def test_original_width_support_and_tails_remain(self):
        self.assertIn('C_setup^-20',protocol()['fixed']);self.assertIn('parent additive ledger retained',protocol()['tails'])
    def test_frozen_record_and_bindings(self):
        old=json.loads((HERE/'evidence.json').read_text());validate(old);compare_records(self.record,old)
        self.assertEqual(len(old['source_hashes']),24);self.assertGreaterEqual(old['inherited_provenance']['records'],80)
    def test_scope_promotions_rejected(self):
        for key,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):validate(r)
    def test_erased_history_and_limits_rejected(self):
        for key in ('inherited_provenance','remaining_obligations'):
            r=deepcopy(self.record)
            if key=='inherited_provenance':r[key]['bindings'].pop(next(iter(r[key]['bindings'])))
            else:r[key]=[]
            with self.assertRaises(ValueError):validate(r)
    def test_exact_seed_and_counterexample_corruption_rejected(self):
        r=deepcopy(self.record);r['reviews']['fresh_canonical_prefix']['orders']['final_regular']='2'
        with self.assertRaises(ValueError):validate(r)
        r=deepcopy(self.record);r['reviews']['arbitrary_diagonal_counterexample']['rows'][0]['log2_diagonal_residual']='0'
        with self.assertRaises(ValueError):validate(r)
    def test_protocol_and_failed_gate_rejected(self):
        r=deepcopy(self.record);r['protocol']['family']='same old fields'
        with self.assertRaises(ValueError):validate(r)
        r=deepcopy(self.record);r['gates']['whole_half_line_operator_bridge']=False
        with self.assertRaises(ValueError):validate(r)
    def test_frozen_output_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
    def test_json_roundtrip(self):validate(json.loads(json.dumps(self.record)))
if __name__=='__main__':unittest.main()
