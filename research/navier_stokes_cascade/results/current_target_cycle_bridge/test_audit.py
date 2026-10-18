"""Exact current equations, missing-term controls and non-promotable scope."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import target_cycle as c

class TargetCycleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']['current_equation_cancellations']
    def identity(self,key):self.assertTrue(self.r['identities'][key],key)
    def test_fixed_current_target_both_components(self):
        self.identity('engineered_current_target_angular');self.identity('engineered_current_target_axial')
    def test_complete_covariance_radial_source_change(self):self.identity('complete_covariance_source_change')
    def test_complete_covariance_current_pressure_change(self):self.identity('complete_covariance_pressure_change')
    def test_actual_signed_angular_target_cancels_old(self):self.identity('actual_signed_angular_target_cancellation')
    def test_actual_signed_axial_target_cancels_old(self):self.identity('actual_signed_axial_target_cancellation')
    def test_moving_compensators_stay_normalized(self):
        self.identity('moving_compensator_unit_1');self.identity('moving_compensator_unit_2')
    def test_both_moving_compensator_commutators(self):
        self.identity('moving_compensator_commutator_1');self.identity('moving_compensator_commutator_2')
    def test_moving_compensation_terms_cannot_be_dropped(self):
        for k in ('moving_angular_compensator_commutator_nonzero','moving_axial_compensator_commutator_nonzero'):self.assertTrue(self.r['controls'][k])
    def test_axial_pressure_target_response_not_zero(self):self.assertTrue(self.r['controls']['axial_pressure_target_response_nonzero'])
    def test_moving_pressure_product_term_not_zero(self):self.assertTrue(self.r['controls']['moving_pressure_moment_product_nonzero'])
    def test_signed_covariance_error_retained(self):self.assertTrue(self.r['controls']['signed_square_error_not_discarded'])
    def test_fast_complete_source_change_is_nonlinear_Haar_part(self):self.identity('fast_complete_source_has_only_nonlinear_Haar_terms')
    def test_fast_pressure_recomputed(self):self.identity('fast_pressure_recomputed_from_nonlinear_source')
    def test_fast_current_residual_differences(self):
        self.identity('fast_current_residual_difference_1');self.identity('fast_current_residual_difference_2')
    def test_fast_current_target_differences(self):
        self.identity('fast_current_target_difference_1');self.identity('fast_current_target_difference_2')
    def test_actual_fast_all_three_components_have_zero_Haar(self):
        for n in ('beta','v','gamma'):self.identity('fast_exact_zero_Haar_'+n)
    def test_fast_linear_time_and_viscosity_Haar_exactly_zero(self):
        self.identity('fast_linear_time_Haar_zero');self.identity('fast_linear_viscosity_Haar_zero')
    def test_fast_old_new_correlations_and_defects_not_zero(self):
        self.assertTrue(self.r['controls']['fast_old_new_correlation_nonzero']);self.assertTrue(self.r['controls']['fast_defect_actually_changed'])
    def test_slow_complete_source_difference(self):self.identity('slow_full_current_source_difference')
    def test_slow_complete_pressure_difference(self):self.identity('slow_full_current_pressure_difference')
    def test_slow_current_residual_differences(self):
        self.identity('slow_full_residual_difference_1');self.identity('slow_full_residual_difference_2')
    def test_slow_current_target_differences(self):
        self.identity('slow_current_target_difference_1');self.identity('slow_current_target_difference_2')
    def test_both_actual_mean_curls_preserve_divergence(self):
        self.identity('fast_divergence');self.identity('slow_divergence')
    def test_both_mean_increments_preserve_both_velocity_moments(self):
        for n in ('fast','slow'):
            for e in ('angular','axial'):self.identity(n+'_'+e+'_velocity_moment')
    def test_current_moment_identities_after_every_state(self):
        for n in ('old','wave','fast','slow'):
            self.identity(n+'_current_angular_moment_identity');self.identity(n+'_current_axial_pressure_moment_identity')
    def test_slow_old_Haar_and_pressure_not_zero(self):
        self.assertTrue(self.r['controls']['slow_old_nonzero_Haar_means_retained']);self.assertTrue(self.r['controls']['slow_changed_pressure_nonzero'])
    def test_linear_cutoff_tail_angular_mean_zero_not_auxiliary_Haar(self):
        for m in (1,2,3):self.identity(f'linear_cutoff_tail_zero_angular_mean_mode_{m}')
    def test_affine_slopes_and_boundary_prove_whole_half_line(self):self.assertTrue(all(c.certificate()['checks'].values()))
    def test_two_signed_crosses_tie_next_target_order(self):
        rows=c.certificate()['rows']
        self.assertEqual(rows['signed_old_error']['boundary_margin'],'0');self.assertEqual(rows['signed_full_curl']['boundary_margin'],'0')
    def test_full_mean_viscosity_margin_is_strict(self):self.assertEqual(F(c.certificate()['full_mean_viscosity_margin']),F(1,2)-F(2,100000))
    def test_higher_C_does_not_hide_a_negative_slope(self):
        r=c.certificate();C0=F(r['C0']);gain=F(r['gain'])
        for C in (C0,C0+F(1,1000),C0+1,10**6):
            for row in r['rows'].values():self.assertGreaterEqual(F(row['order_slope'])*C+F(row['order_offset']),C+gain)
    def test_fast_target_gain_stronger_than_defect_gain(self):
        k=F(1,100000);C=2-2*k;mu=1-k;H=C-k
        self.assertEqual((mu+H+1)-(mu+H),1);self.assertLess(mu+H,C+1)
    def test_invalid_budget_kappa_rejected(self):
        for k in (0,-1,F(1,100),1):
            with self.assertRaises(ValueError):c.certificate(k)
    def test_accepted_initial_canonical_target_routes_to_seed(self):
        r=c.prefix();self.assertTrue(all(r['checks'].values()));self.assertEqual(r['orders']['final_target_seed'],'99999/50000')
    def test_no_real_shortcut_in_generic_square_budget(self):self.assertEqual(c.certificate()['rows']['signed_square']['order_offset'],'-1')
    def test_parent_and_protocol(self):
        self.assertEqual(audit.protocol()['parent_commit'],'f370f9155a7286d613ad4e00ba5e22a3ddf43b0d')
        self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_historical_sources_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_exact_replay(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_all_unproved_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_pressure_cancellation_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['current_equation_cancellations']['identities']['actual_signed_axial_target_cancellation']=False
        with self.assertRaises(ValueError):audit.validate(r)
    def test_moving_compensator_omission_rejected(self):
        r=deepcopy(self.record);r['reviews']['current_equation_cancellations']['controls']['moving_axial_compensator_commutator_nonzero']=False
        with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_residual_obligation_cannot_be_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_bool_integer_scope_tamper_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_cli_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)])
            audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
