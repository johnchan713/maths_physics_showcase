"""Exact factored-target, nonlinear chart, cutoff and same-root regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,profile_fixture,original_fixture,fixed_profile
from bounds import P,ChartBudget,TargetBudget,CutoffBudget,from_original

class PrimaryTargetCutoffTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        vals=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(vals,key);self.assertTrue(all(vals),key)
    def test_smooth_factored_profile_reconstructs_original_target(self):self.control('outer_factor_reconstructs_same');self.control('inner_cancellation_leaves_smooth')
    def test_mixed_quotient_edge_jets_on_whole_collar(self):self.control('all_mixed_outer_collar_quotient_jets_enclosed')
    def test_cubic_pole_is_not_falsely_declared_bounded(self):self.control('quotient_itself_keeps_cubic_outer_pole')
    def test_original_flat_constants_not_replaced(self):self.control('replacing_inner_flat_constant_changes');self.control('original_width_square_is_retained')
    def test_sharper_pole_supplies_original_root_allowance(self):self.control('sharper_pole_fits_original_root')
    def test_exact_implicit_log_coordinate_bounds(self):self.control('implicit_log_coordinate_majorant')
    def test_eta_and_original_band_scale_prefactor_jets(self):self.control('eta_and_scale_factor_majorants')
    def test_full_target_composition_through_order_four(self):self.control('original_band_target_composition_bound')
    def test_full_loghalf_composition_through_order_four(self):self.control('loghalf_composition_bound')
    def test_nonlinear_chart_terms_are_not_dropped(self):self.control('nonlinear_chart_terms_cannot_be_omitted')
    def test_loghalf_zero_order_slot_is_not_claimed_bounded(self):self.control('loghalf_order_zero_is_unused')
    def test_auxiliary_dual_covector_cutoff_jets(self):self.control('dual_auxiliary_covector_cutoff_jets_enclosed')
    def test_original_S_cubed_mesh_scale_retained(self):self.control('original_S_cubed_mesh_jets_enclosed')
    def test_temporal_fast_time_length_cancellation_exact(self):self.control('relative_pulse_cutoff_cancels_fast_time_length')
    def test_complete_cutoff_has_finite_polynomial_scale_degree(self):self.control('complete_cutoff_polynomial_degree')
    def test_slow_partition_is_in_wave_norm_once_and_outside_H(self):self.control('slow_partition_is_included_once');self.assertTrue(self.record['flags']['slow_partition_included_once_in_wave_norm_and_kept_outside_H'])
    def test_fixed_radius_and_width_not_reselected(self):self.control('original_fixed_geometry_radius');self.control('original_width_and_radius_retained')
    def test_original_band_and_tangent_normal_gates_retained(self):self.control('original_tail_and_normal_tangent')
    def test_all_five_root_norm_slots_are_derived(self):self.control('target_and_loghalf_slots_derived');self.control('H_wave_and_complete_curl_slots_derived')
    def test_U_chart_transports_target_loghalf_and_cutoff_wave(self):self.control('U_chart_transports_target_loghalf_and_cutoff_wave')
    def test_value_gaps_and_same_root_domain_retained(self):self.control('all_original_root_value_gaps')
    def test_complete_curl_and_signed_recipes_are_finite(self):self.control('joined_amplitude_curl_and_signed')
    def test_invalid_width_or_radial_hull_rejected(self):
        for args in ((0,1),(-1,1),(1,0),(1,2)):
            with self.assertRaises(ValueError):ChartBudget(*args)
        with self.assertRaises(TypeError):ChartBudget(0.5,1)
    def test_invalid_fixed_profile_seminorms_rejected(self):
        with self.assertRaises(TypeError):TargetBudget(object(),lambda j:1)
        with self.assertRaises(ValueError):TargetBudget(ChartBudget(1,1),lambda j:-1).target(0)
        with self.assertRaises(TypeError):TargetBudget(ChartBudget(1,1),lambda j:0.5).target(0)
    def test_invalid_cutoff_input_and_jet_orders_rejected(self):
        c=ChartBudget(1,1);kw=dict(dyadic=fixed_profile,mesh=fixed_profile,transverse=fixed_profile,temporal=fixed_profile)
        with self.assertRaises(ValueError):CutoffBudget(c,**kw,radius=0)
        with self.assertRaises(TypeError):CutoffBudget(c,**kw).normalized(0,1)
        for n in (-1,Q(1,2)):
            with self.assertRaises((ValueError,TypeError)):profile_fixture()[-1].target(n)
        with self.assertRaises(ValueError):c.increment('wrong',1)
    def test_unknown_band_is_not_reselected(self):
        p,sp,ell,kw=original_fixture()
        with self.assertRaises(ValueError):from_original(p,sp,16,**kw)
    def test_all_controls_and_historical_sources_match(self):self.assertTrue(all(self.r['checks'].values()));self.assertTrue(all(self.record['gates'].values()))
    def test_eight_frozen_dependencies_and_current_source_hashes(self):self.assertEqual(len(self.record['dependency_evidence_sha256']),8);self.assertEqual(self.record['source_hashes'],audit.sources())
    def test_no_numerical_giant_profile_or_full_residual_promotion(self):
        for k in ('actual_giant_profile_seminorms_full_background_or_value_gate_numerically_certified','unconditional_pulse_full_physical_residual_infinite_field_or_blowup_proved'):
            self.assertIs(self.record['flags'][k],False)
        self.assertTrue(self.record['flags']['all_pressure_memory_axial_means_reserved_supports_and_residual_terms_retained'])
    def test_tampered_scope_recipe_or_provenance_rejected(self):
        for field in ('flags','review','source_hashes'):
            r=deepcopy(self.record);r[field]['tampered']=True
            with self.assertRaises(ValueError):audit.validate(r)
    def test_frozen_evidence_reproduces(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()))
    def test_cli_uses_fresh_output_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/'replay.json';audit.main(['--verify-record',str(Path(__file__).parent/'evidence.json'),'--output',str(dest)])
            self.assertEqual(json.loads(dest.read_text()),self.record)
            with self.assertRaises(SystemExit):audit.main(['--output',str(dest)])
if __name__=='__main__':unittest.main()
