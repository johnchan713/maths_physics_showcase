"""Exact current residual/pressure controls and strict conditional-scope guards."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import residual_ledger as r

class CurrentResidualTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.x=cls.record['reviews']['exact_complete_residual_ledger']
    def identity(self,k):self.assertTrue(self.x['identities'][k],k)
    def test_base_divergence_free(self):self.identity('base_divergence_free')
    def test_old_full_cutoff_commutator(self):self.identity('old_complete_cutoff_commutator')
    def test_current_regular_source_modes_are_the_inverted_modes(self):self.identity('current_source_covers_each_inverted_mode')
    def test_old_complete_residual_split(self):self.identity('old_residual_split_is_exact')
    def test_full_current_regular_plus_new_tail_ledger(self):self.identity('complete_current_sourced_residual_ledger')
    def test_both_full_viscous_cutoff_commutators(self):
        for m in (1,2):self.identity(f'full_cutoff_commutator_mode_{m}')
    def test_cut_potential_product_rule_with_full_curl(self):
        for m in (1,2):self.identity(f'cut_potential_equals_product_plus_curl_mode_{m}')
    def test_complete_cut_velocity_divergence_zero(self):
        for m in (1,2):self.identity(f'complete_cut_divergence_zero_mode_{m}')
    def test_cutoff_curl_cancels_naive_cut_divergence(self):
        for m in (1,2):self.identity(f'cutoff_curl_divergence_cancels_mode_{m}')
    def test_naive_cut_velocity_is_not_divergence_free(self):self.assertTrue(self.x['nonzero_omission_controls']['naive_cut_velocity_is_not_divergence_free'])
    def test_every_cutoff_commutator_group_nonzero(self):
        for k in ('temporal','mixed_viscosity','cutoff_laplacian','pressure','cutoff_curl_operator'):self.assertTrue(self.x['nonzero_omission_controls'][k],k)
    def test_both_old_new_crosses_and_new_square_nonzero(self):
        for k in ('both_old_new_cross_orientations_nonzero','new_wave_square_nonzero'):self.assertTrue(self.x['nonzero_omission_controls'][k])
    def test_old_and_uncancelled_source_tails_retained_nonzero(self):
        for k in ('old_tail_retained_nonzero','uncancelled_source_cutoff_retained_nonzero'):self.assertTrue(self.x['nonzero_omission_controls'][k])
    def test_recomputed_mean_pressure_preserves_nonzero_residual(self):self.identity('recomputed_mean_pressure_does_not_change_nonzero_residual')
    def test_full_mean_pressure_change_is_nonzero_gradient(self):
        self.identity('full_mean_pressure_difference_is_gradient');self.assertTrue(self.x['nonzero_omission_controls']['actual_mean_pressure_gradient_change_nonzero'])
    def test_Pi2_axial_memory_is_not_an_erasable_gauge(self):
        self.identity('same_full_base_pressure_Pi2_cancels_only_in_difference');self.assertTrue(self.x['nonzero_omission_controls']['Pi2_axial_gradient_nonzero'])
    def test_actual_mean_increment_divergence_zero(self):self.identity('actual_mean_increment_divergence_zero')
    def test_actual_mean_wave_residual_increment(self):self.identity('actual_mean_update_nonzero_residual_exact')
    def test_mean_update_keeps_nonzero_tail_ledger(self):self.identity('mean_update_keeps_exact_nonzero_tail_ledger')
    def test_actual_mean_wave_increment_nonzero(self):self.assertTrue(self.x['nonzero_omission_controls']['actual_mean_wave_change_nonzero'])
    def test_pressure_only_cutoff_commutator(self):
        self.identity('pressure_only_cutoff_commutator');self.assertTrue(self.x['nonzero_omission_controls']['pressure_only_cutoff_term_nonzero'])
    def test_actual_temporal_factor_has_no_spatial_curl(self):self.identity('temporal_cutoff_has_no_spatial_curl')
    def test_actual_temporal_factor_has_no_pressure_spatial_commutator(self):self.identity('temporal_cutoff_has_no_spatial_pressure_commutator')
    def test_temporal_tail_is_only_time_product(self):
        self.identity('temporal_cutoff_tail_is_only_time_product');self.assertTrue(self.x['nonzero_omission_controls']['temporal_time_product_nonzero'])
    def test_generic_spatial_peak_cannot_be_declared_Gaussian_flat(self):
        q=r.gaussian_review()['spatial_peak_not_Gaussian_tail'];self.assertTrue(q['nonzero_spatial_derivative_at_temporal_peak']);self.assertGreater(F(q['log_growth_slope']),0)
    def test_conjugate_reality_residual_and_tails(self):
        self.identity('conjugate_complete_residual_reality');self.identity('conjugate_exact_tail_reality')
    def test_highest_generated_harmonic_four_is_retained(self):
        self.assertTrue(self.x['nonzero_omission_controls']['highest_generated_harmonic_four_present']);self.assertIn(4,self.x['new_harmonics'])
    def test_all_six_complete_cut_covariance_entries(self):
        for i in range(3):
            for j in range(i,3):self.identity(f'complete_cut_covariance_{i}{j}')
    def test_whole_half_line_residual_budget(self):self.assertTrue(all(r.cycle_certificate()['checks'].values()))
    def test_seven_limiting_groups_tie(self):
        rows=r.cycle_certificate()['rows'];self.assertEqual(sum(v['margin_slope']=='0' and v['boundary_margin']=='0' for v in rows.values()),7)
    def test_square_and_slow_groups_have_strict_margins(self):
        rows=r.cycle_certificate()['rows']
        for k in ('sourced_complete_square','signed_complete_square','actual_slow_mean_wave'):self.assertGreater(F(rows[k]['boundary_margin']),0)
    def test_affine_certificate_checks_large_B_without_new_domain(self):
        c=r.cycle_certificate()
        for B in (F(c['B0']),10,10**12):
            for v in c['rows'].values():self.assertGreaterEqual(F(v['order_slope'])*B+F(v['order_offset']),B+F(c['gain']))
    def test_flat_certificate_zero_loss_and_zero_polynomial(self):
        c=r.flat_certificate(0,0,0,1);self.assertEqual(c['prefactor'],'1');self.assertEqual(c['quadratic_log_constant'],'0')
    def test_flat_certificate_retains_small_Gaussian_constant_cost(self):
        self.assertGreater(F(r.flat_certificate(1,1,4,F(1,1000))['prefactor']),F(r.flat_certificate(1,1,4,1)['prefactor']))
    def test_fractional_logarithmic_degree_rounded_up(self):self.assertEqual(r.flat_certificate(1,2,F(7,2),1)['series_order'],2)
    def test_exact_completed_square_controls(self):self.assertTrue(all(x['exact_completed_square_controls'] for x in r.gaussian_review()['rows']))
    def test_flat_certificate_invalid_inputs_rejected(self):
        for a in ((-1,0,0,1),(0,-1,0,1),(0,0,-1,1),(0,0,0,0)):
            with self.assertRaises(ValueError):r.flat_certificate(*a)
    def test_flat_certificate_float_and_bool_not_exact_inputs(self):
        for a in ((1.,0,0,1),(0,0,False,1)):
            with self.assertRaises(TypeError):r.flat_certificate(*a)
    def test_normalized_source_gain_is_not_physical_flatness(self):
        rows=r.gaussian_review()['normalized_source_obstruction'];self.assertTrue(all(F(x['log2_normalized_ratio'])>0 for x in rows[1:]))
    def test_parent_and_protocol(self):
        self.assertEqual(audit.protocol()['parent_commit'],'62cce31db067535582f9b22d4ea0727a3c158f04');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_historical_source_bytes_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_exact_replay(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_all_unproved_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                x=deepcopy(self.record);x['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(x)
    def test_pressure_only_commutator_tamper_rejected(self):
        x=deepcopy(self.record);x['reviews']['exact_complete_residual_ledger']['identities']['pressure_only_cutoff_commutator']=False
        with self.assertRaises(ValueError):audit.validate(x)
    def test_required_initial_input_obligation_cannot_be_removed(self):
        x=deepcopy(self.record);x['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(x)
    def test_bool_integer_scope_tamper_rejected(self):
        x=deepcopy(self.record);x['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(x)
    def test_fresh_cli_replay_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)])
            audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
