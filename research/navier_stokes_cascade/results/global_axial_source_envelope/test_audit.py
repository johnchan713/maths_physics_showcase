"""Mixed-order, global-source, full viscosity and scope controls."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import envelopes as w

class SourceEnvelopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def mixed(self,k):self.assertTrue(self.r['mixed_and_reference']['checks'][k])
    def symbolic(self,k):self.assertTrue(self.r['symbolic']['checks'][k])
    def test_weighted_pullback_all_twenty_terms(self):self.mixed('weighted_pullback_all_twenty_terms_bounded_by_H_over_N');self.assertEqual(len(self.r['mixed_and_reference']['pullback_weight_rows']),20)
    def test_small_weighted_not_small_radial_norm(self):self.mixed('H_over_N_less_than_one');self.assertTrue(self.r['pullback']['checks']['third_pure_phase_loss_is_N_squared'])
    def test_exact_all_thirty_raw_chain_controls(self):self.assertEqual(len(self.r['pullback']['rows']),30);self.assertTrue(self.r['pullback']['checks']['all_thirty_independent_chain_controls'])
    def test_same_repair_mixed_J3_below_one(self):self.mixed('final_repair_C3_mixed_is_below_one')
    def test_same_N_native_J3_envelope(self):self.mixed('original_and_repair_fit_same_mixed_bound')
    def test_H70_not_H64_relabel(self):self.mixed('fixed_frequency_upper_gives_H70')
    def test_first_selective_log_jet_no_N(self):self.mixed('selective_first_log_jet_does_not_need_N')
    def test_second_selective_log_jet_keeps_N(self):self.mixed('selective_second_log_jet_retains_N')
    def test_second_log_difference_keeps_lower_order(self):self.mixed('log_second_difference_includes_lower_derivative')
    def test_actual_unit_mass_repair_width(self):self.mixed('actual_unit_mass_width_point_one_third_jet_cap')
    def test_derivative_cap_absorbed_at_actual_C_floor(self):self.mixed('bump_derivative_constant_absorbed_by_C12')
    def test_reference_shape_three_piece_caps(self):self.assertEqual(self.r['mixed_and_reference']['reference_R0_derivative_caps'],['11','353','13440','3040000']);self.mixed('R0_three_joined_pieces_fit_one_e7')
    def test_reference_pulse_end_bumps_not_unit_mass_repair_bumps(self):self.mixed('actual_width_point_three_reference_bump_fits_one_e103')
    def test_complete_reference_axial_pulse(self):self.mixed('full_reference_axial_pulse_fits_one_e110')
    def test_reference_taper_jets(self):self.mixed('reference_taper_three_derivatives_fit_one_e7')
    def test_actual_H_dominates_later_reference(self):self.mixed('actual_H_above_two_to_1024_dominates_reference')
    def test_exact_VX_identity(self):self.symbolic('exact_VX_identity')
    def test_exact_VXX_lowers_derivative_order(self):self.symbolic('exact_VXX_needs_only_U_third_total_jets')
    def test_full_shifted_axial_operator(self):self.symbolic('full_shifted_axial_operator')
    def test_omitted_axial_coefficient_derivatives_detected(self):self.symbolic('omitted_axial_coefficient_derivatives_detected')
    def test_log_second_chain(self):self.symbolic('log_second_retains_lower_derivative')
    def test_reserved_plateau_pressure_source_nonzero(self):self.symbolic('reserved_plateau_keeps_nonzero_pressure_source')
    def test_reciprocal_closed_h_interval(self):self.assertTrue(self.r['axial_coefficients']['checks']['reciprocal_L_C2_strictly_below_two'])
    def test_shifted_operator_all_three_coefficients(self):self.assertTrue(all(self.r['axial_coefficients']['checks'].values()))
    def test_all_ninety_independent_source_controls(self):self.assertEqual(self.r['manufactured']['case_count'],90);self.assertTrue(self.r['manufactured']['checks']['all_ninety_exact_pointwise_controls'])
    def test_source_value_and_eta_jets_checked(self):self.assertTrue(all(r['checks']['exact_pointwise_source_and_first_angular_derivative'] for r in self.r['manufactured']['rows']))
    def test_axis_and_both_angular_endpoints(self):self.assertTrue(self.r['manufactured']['checks']['both_axis_and_angular_endpoints_retained'])
    def test_pointwise_viscosity_nonzero(self):self.assertTrue(self.r['manufactured']['checks']['viscosity_is_not_a_pointwise_zero_term'])
    def test_viscosity_bound_slope_retained(self):
        a=w.source_bounds(F(1,200),0,3,5,7,11);b=w.source_bounds(F(1,200),F(7,3),3,5,7,11)
        self.assertEqual(b['w0_C1']-a['w0_C1'],F(14,3)*a['VXX_C1'])
    def test_absolute_mass_is_full_Xv(self):
        r=w.source_bounds(F(1,200),1,11,5,7,13);self.assertEqual(r['absolute_source_mass'],11*r['w0_C1']);self.assertEqual(r['absolute_source_first_mass'],F(121,2)*r['w0_C1'])
    def test_global_KU_uses_later_pulse(self):self.assertEqual(w.global_bounds(F(1,200),1,11,1,7,13,101,303)['U_C3'],352)
    def test_global_KU_uses_axis(self):self.assertEqual(w.global_bounds(F(1,200),1,11,1,1000,13,101,303)['U_C3'],1000)
    def test_global_KU_uses_modulation(self):self.assertEqual(w.global_bounds(F(1,200),1,11,1,7,1000,101,303)['U_C3'],2000)
    def test_axis_radial_envelope_not_divided_by_X(self):
        r=w.global_bounds(F(1,200),1,11,1,1000,13,1,1);self.assertEqual(r['U_X_C2'],1000);self.assertEqual(r['U_XX_C1'],1000)
    def test_log_to_X_inverse_radius_costs(self):
        r=w.global_bounds(F(1,200),1,11,F(1,2),0,13,7,11);self.assertEqual(r['U_X_C2'],14);self.assertEqual(r['U_XX_C1'],44)
    def test_global_v_angular_norm_retains_factor_36(self):self.assertEqual(F(self.r['auxiliary_global_envelope']['v_C2']),36*352)
    def test_axial_source_formula_coefficients(self):self.assertEqual(w.scalar_axial_bound(2,3,5,7),64*4*7+704*2*5+1664*3)
    def test_zero_U_has_zero_source(self):self.assertEqual(w.source_bounds(F(1,200),1,11,0,0,0)['w0_C1'],0)
    def test_invalid_source_inputs(self):
        for args in ((0,1,1,1,1,1),(F(1,50),1,1,1,1,1),(F(1,200),-1,1,1,1,1),(F(1,200),1,0,1,1,1),(F(1,200),1,1,-1,1,1)):
            with self.assertRaises(ValueError):w.source_bounds(*args)
    def test_invalid_global_geometry(self):
        with self.assertRaises(ValueError):w.global_bounds(F(1,200),1,1,2,7,13,7,11)
    def test_exact_rational_inputs_only(self):
        for q in (True,.1):
            with self.assertRaises(TypeError):w.source_bounds(F(1,200),1,11,q,7,13)
    def test_parent_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'dfcf62462a6a10da6313be53638f3ab2b25d9583');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_historical_source_bytes(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_record(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_no_scope_promotion(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_scope_type_mutation(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_swirl_obligation_not_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_pointwise_budget_tamper(self):
        r=deepcopy(self.record);r['reviews']['auxiliary_global_envelope']['w0_C1']='0'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
