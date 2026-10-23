"""Full compact mean-source quotient and pressure/viscous regression tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from bounds import rates,reduced_envelope,A1_request,certificate
from controls import review

class CompactQuotientTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_full_radial_vector_viscosity_factor(self):self.control('complete_radial_vector_viscosity_is_DrD1')
    def test_full_axial_scalar_viscosity_factor(self):self.control('complete_axial_scalar_viscosity_is_D1Dr')
    def test_complete_gr_reduction_with_K0(self):self.control('full_gr_is_reduced_r0_plus_Dr_K0')
    def test_full_Ez_difference_keeps_all_time_pressure_and_viscous_parts(self):self.control('full_Ez_retains_time_pressure_flux_and_both_viscosities')
    def test_entire_current_Ez_has_explicit_complete_D1_quotient_flux(self):self.control('entire_current_Ez_minus_Dz_H1_is_one_complete_D1_flux')
    def test_same_actual_potential_gives_exact_incompressibility(self):self.control('actual_complete_potential_is_exactly_divergence_free')
    def test_time_of_gamma_is_in_full_D1_range(self):self.control('time_of_actual_gamma_stays_in_D1_range')
    def test_both_viscosities_of_gamma_are_in_full_D1_range(self):self.control('both_axial_viscosities_of_actual_gamma_stay_in_D1_range')
    def test_radial_viscous_pressure_memory_is_not_deleted_from_field(self):
        self.control('K0_pressure_memory_is_not_zero');self.control('radial_viscous_pressure_memory_annihilated_only_after_A1Dz')
    def test_exact_weighted_characteristic_mismatch_and_boundary_formula(self):self.control('exact_full_drift_weight_and_boundary_rule')
    def test_A_e_annihilates_only_the_matching_radial_derivative(self):self.control('A_e_annihilates_compact_D_e_potential')
    def test_pressure_primitive_restores_compact_potential(self):self.control('T_e_D_e_restores_pressure_potential_not_zero')
    def test_noncompact_boundary_term_blocks_local_annihilation(self):self.control('noncompact_boundary_cannot_be_erased')
    def test_same_cutoff_gives_projection_and_TA_zero(self):
        self.control('same_cutoff_endpoint_mass_gives_A_e_projection');self.control('T_e_A_e_zero_uses_same_chi_and_not_a_second_cutoff')
    def test_second_cutoff_mutation_detected(self):self.control('different_cutoff_in_same_operator_identity_is_detected')
    def test_full_actual_A1_velocity_is_essential(self):
        self.control('actual_gd_minus_A1_gd_is_in_kernel_A1');self.control('dropping_actual_A1_velocity_breaks_kernel_identity')
    def test_moving_cutoff_both_slow_commutators_and_time_sign(self):
        self.control('moving_chi_Z_commutator_is_retained');self.control('moving_chi_T_commutator_is_retained');self.control('tstar_commutator_has_negative_epsilon_sign')
    def test_free_Dz_commutation_is_rejected(self):self.control('commuting_Dz_past_A1_without_chi_term_changes_result')
    def test_minus_S_pressure_memory_changes_axial_flux(self):self.control('dropping_minus_S_local_pressure_memory_changes_A1_axial_flux')
    def test_complete_current_compensated_pressure_and_moving_moment(self):
        self.control('complete_current_pressure_is_K0_plus_compensated_primitive');self.control('radial_compensation_not_replaced_by_pointwise_source_zero');self.control('whole_moving_pressure_moment_derivative_keeps_c_rho_jet')
    def test_comparison_swirl_respects_divergence_and_both_velocity_moments(self):
        self.control('pure_swirl_is_divergence_free_beta_gamma_zero');self.control('pure_swirl_exact_Haar_velocity_moments_both_zero')
    def test_surviving_centripetal_fundamental_and_square_modes(self):
        self.control('full_centripetal_fundamental_is_exact_resonant_source');self.control('centripetal_square_cannot_cancel_fundamental');self.control('actual_parent_frequency_is_nonzero')
    def test_swirl_pressure_defect_kept_and_compensation_cannot_erase_fundamental(self):
        self.control('Haar_radial_pressure_change_is_nonzero_and_must_be_compensated');self.control('same_scalar_pressure_compensator_changes_no_fundamental')
    def test_full_positive_polynomial_and_raw_epsilon_budget(self):
        self.control('full_positive_polynomial_and_exact_rational_envelope_agree');self.control('independent_full_pressure_flux_B_and_A1_B_plus_one_budget')
    def test_explicit_Dz_epsilon_Qcirc_and_n_plus_three_input_cost(self):
        self.control('A1_Qcirc_factor_two_and_explicit_Dz_epsilon_retained');self.control('whole_pressure_flux_derivative_has_n_plus_three_input_budget')
    def test_covariance_and_pressure_norm_mutations_change_budget(self):
        self.control('complete_covariance_rr_rz_tt_zz_not_dropped');self.control('full_current_pressure_operator_changes_A1_budget')
    def test_conditional_A1_gain_two_kappa_not_full_row_promotion(self):self.control('r0_remains_B_and_A1_only_gains_two_kappa')
    def test_normalized_A1_direct_operator_with_zero_integrations(self):self.control('normalized_A1_uses_direct_operator_no_oscillatory_integration')
    def test_physical_A1_H1_extra_derivative_and_strict_margin(self):
        self.control('physical_raw_H1_keeps_extra_Dz_and_strict_margin');self.control('illustrative_physical_A1_p30000002_H1_J60000008')
    def test_distinct_stages_domains_and_older_A0_cost(self):
        self.control('postwave_A1_and_final_A0_fields_not_mixed');self.control('original_6q_KA1_and_12q_final_radial_domains_retained');self.control('old_physical_A0_cost_not_erased')
    def test_regular_E_and_actual_input_gate_remain(self):self.control('whole_regular_E_and_actual_input_gate_unchanged')
    def test_invalid_inputs_rejected(self):
        for k in (0,Q(1,4)):
            with self.assertRaises(ValueError):rates(k)
        with self.assertRaises(ValueError):reduced_envelope(0,{})
        with self.assertRaises(ValueError):A1_request(0,physical_power=-1)
        with self.assertRaises(ValueError):certificate(-1)
    def test_review_cache_cannot_alias_frozen_evidence(self):
        r=review();r['checks'].clear();self.assertTrue(review()['checks'])
    def test_six_history_records_and_frozen_evidence_reproduce(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),6)
        self.assertTrue(all(self.record['gates'].values()))
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_protocol_provenance_and_remaining_gate_mutations_rejected(self):
        for field in ('flags','protocol','dependency_evidence_sha256','remaining_obligations'):
            r=deepcopy(self.record)
            if field=='flags':r[field]['whole_Ez_pressure_or_complete_regular_E_order_promoted']=True
            elif field=='protocol':r[field]['pressure']='discard K0'
            elif field=='dependency_evidence_sha256':r[field][next(iter(r[field]))]='wrong'
            else:r[field].pop()
            with self.assertRaises(ValueError):audit.validate(r)
    def test_frozen_outputs_never_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')
if __name__=='__main__':unittest.main()
