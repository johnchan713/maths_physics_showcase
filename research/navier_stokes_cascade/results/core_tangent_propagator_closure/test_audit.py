#!/usr/bin/env python3
"""Whole-vector energy, endpoint frames, canonical zero entry and inverse join."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,fixture,psd2
from budgets import normal,transfer,certificate,tangent_polynomial,forced_from_frame,residence_upper


class TangentClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        rows=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(rows);self.assertTrue(all(rows))

    def test_energy_gap_certifies_all_complex_vectors(self):self.control('whole_complex_vector_logarithmic_energy_bound')
    def test_omitting_damping_deviation_fails(self):self.control('omitting_damping_deviation_breaks_energy_bound')
    def test_exact_frame_Gram_and_singular_values(self):self.control('exact_Gram_and_singular_value_squares')
    def test_specified_left_inverse_is_used_on_actual_tangent_plane(self):self.control('specified_left_inverse_restricts_to_identity')
    def test_endpoint_frame_conversion_cannot_be_omitted(self):self.control('coordinate_identity_can_expand_physical_tangent_norm')
    def test_full_two_time_gauge_preserves_physical_map(self):self.control('two_time_physical_tangent_maps_are_identical_under_full_gauge')
    def test_zero_entrance_is_exactly_invariant(self):self.control('zero_entrance_stays_zero_in_either_frame')
    def test_homogeneous_primary_mismatch_is_preserved(self):self.control('resetting_growing_coordinate_changes_initial_physical_pulse')
    def test_complete_time_dependent_generator_and_forcing_transform(self):
        self.control('exact_time_dependent_generator_conjugacy');self.control('forcing_coordinates_transform_exactly')

    def test_scalar_damping_does_not_require_commuting_coefficients(self):
        self.control('coefficient_evolutions_need_not_commute');self.control('scalar_damping_factors_through_noncommuting_evolution')

    def test_all_nonzero_harmonics_keep_extra_contraction(self):self.control('extra_actual_damping_never_amplifies')
    def test_normal_floor_and_mode_residence(self):
        self.control('normal_floor_supplies_positive_damping_floor');self.control('residence_gain_for_high_modes_not_fundamental_modes')

    def test_original_error_constants_are_used_without_enlarging_length(self):self.control('original_CE_plus_Cd_is_the_entire_integrated_cost')
    def test_exact_constant_recipe_is_scale_independent(self):
        self.control('exact_exponential_upper_recipe_without_float');self.control('tangent_bound_is_constant_in_S')

    def test_accepted_geometric_length_replaces_independent_input(self):self.control('original_L_le_S_closes_whole_path_length_input')
    def test_normal_denominator_and_source_callback_unchanged(self):
        self.control('original_normal_denominator_unchanged');self.control('source_callback_not_estimated_or_changed')

    def test_joined_inverse_retains_all_full_normal_jets(self):self.control('full_normal_and_moving_normal_callbacks_retained')
    def test_joined_mixed_jet_velocity_and_pressure_budgets_are_finite(self):self.control('joined_velocity_and_pressure_are_finite_positive_polynomials')
    def test_no_new_acceptance_or_scope_promotion(self):
        self.control('all_original_gates_and_no_new_acceptance');self.control('no_actual_enclosure_primary_identity_or_residual_promotion')

    def test_join_does_not_reselect_band_or_assume_a_tangent_input(self):
        p,sp,ell,kw=fixture();old=normal.selector.select_tail
        def forbidden(*a,**k):raise AssertionError('Band reselection is forbidden')
        normal.selector.select_tail=forbidden
        try:
            b,f,c=forced_from_frame(p,sp,ell,**kw)
            self.assertTrue(c['tangent_and_length_inputs_derived_from_retained_frame_bounds'])
        finally:normal.selector.select_tail=old
        with self.assertRaises(TypeError):forced_from_frame(p,sp,ell,tangent=transfer.P(1),**kw)

    def test_invalid_old_band_length_or_enlarged_primitive_is_rejected(self):
        p,sp,ell,kw=fixture()
        for pp,ss,ee in ((p,sp,16),(p,dict(sp,length=1),ell),(dict(p,Mz=10**30),sp,ell)):
            with self.assertRaises(ValueError):certificate(pp,ss,ee)

    def test_exponential_recipe_cannot_be_rounded_down(self):
        p,sp,ell,_=fixture();c=certificate(p,sp,ell)
        c['integer_exponential_upper_power']-=1
        with self.assertRaises(ValueError):tangent_polynomial(c)

    def test_invalid_harmonic_and_residence_length_rejected(self):
        p,_,_,_=fixture()
        for length,m in ((1,0),(0,1),(1,1.5),(1,True)):
            with self.assertRaises(ValueError):residence_upper(p,length,m)

    def test_PSD_control_rejects_negative_minor_and_nonsymmetry(self):
        self.assertFalse(psd2([[1,2],[2,1]]));self.assertFalse(psd2([[1,1],[0,1]]))

    def test_four_frozen_records_no_new_analytical_assumptions(self):
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertEqual(self.record['protocol']['new_acceptance_inequalities'],[])
        self.assertEqual(len(self.record['dependency_evidence_sha256']),4);self.assertTrue(all(self.record['gates'].values()))

    def test_new_frozen_record_reproduces_exactly(self):self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)
    def test_scope_and_remaining_obligations_cannot_be_promoted(self):
        for k,v in self.record['flags'].items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)

    def test_frozen_output_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
